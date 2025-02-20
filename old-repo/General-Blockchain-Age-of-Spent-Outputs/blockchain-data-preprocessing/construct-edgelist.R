# install.packages("data.table")
# install.packages("RSQLite")
# install.packages("DBI")

library(data.table)
library(RSQLite)
library(DBI)

data.dir <- ""
# Input data directory here, with trailing "/"

source("https://gist.githubusercontent.com/jeffwong/5925000/raw/bf02ed0dd2963169a91664be02fb18e45c4d1e20/sqlitewritetable.R")
# From https://gist.github.com/jeffwong/5925000
# Modifies RSQLite's sqliteWriteTable function so as to reject duplicates

con <- DBI::dbConnect(RSQLite::SQLite(), paste0(data.dir, "tx-graph-node-indices.db"))

DBI::dbExecute(con, "PRAGMA max_page_count = 4294967292;")
# Allows SQL database files up to 4 TB. See:
# https://stackoverflow.com/questions/16685016/sqlite3-operationalerror-database-or-disk-is-full-on-lustre

# DBI::dbExecute(con, "PRAGMA temp_store = 2;")
# Put temp file in RAM:
# https://stackoverflow.com/a/19259699

DBI::dbExecute(con, "CREATE TABLE nodes (
node TEXT,
node_index INTEGER PRIMARY KEY AUTOINCREMENT,
unique(node)
)")


DBI::dbWriteTable(con, "edgelist", 
  data.frame(origin = character(0), destination = character(0), value = numeric(0),
    block_height = integer(0), stringsAsFactors = FALSE))

tx.graph.files <- list.files(paste0(data.dir, "tx_graphs/"))
tx.graph.files <- tx.graph.files[grepl("^tx_graph.+rds$", tx.graph.files)]
tx.graph.files <- sort(tx.graph.files)

tx.graph.indexed <- vector("list", length(tx.graph.files))
names(tx.graph.indexed) <- tx.graph.files

for (file.iter in tx.graph.files) {
  
  tx.graph.chunk <- readRDS(paste0(data.dir, "tx_graphs/", file.iter))
  
  tx.graph.chunk <-
    rbind(
      data.table(origin = paste0(tx.graph.chunk$incoming$origin.txid, "-", 
        formatC(tx.graph.chunk$incoming$origin.position, width = 4, format = "f", flag = "0", digits = 0)),
        destination = tx.graph.chunk$incoming$txid, 
        value = NA_real_,
        block_height = as.integer(tx.graph.chunk$incoming$block_height), stringsAsFactors = FALSE),
      data.table(origin = tx.graph.chunk$outgoing$txid,
        destination = paste0(tx.graph.chunk$outgoing$txid, "-", 
          formatC(tx.graph.chunk$outgoing$position, width = 4, format = "f", flag = "0", digits = 0)),
        value = tx.graph.chunk$outgoing$value,
        block_height = as.integer(tx.graph.chunk$outgoing$block_height), stringsAsFactors = FALSE)
    )
  
  DBI::dbWriteTable(con, "edgelist", 
    tx.graph.chunk, append = TRUE)
  
  cat(file.iter, base::date(), "\n")
  
  if (nrow(tx.graph.chunk) == 0) {next}
  
  new.nodes <- unique(c(tx.graph.chunk$origin, tx.graph.chunk$destination))
  
  nodes.to.insert <- data.frame(node = new.nodes, node_index = NA, stringsAsFactors = FALSE)
  
  mysqliteWriteTable(con, "nodes", 
    nodes.to.insert, append = TRUE, row.names = FALSE, ignore = TRUE)
  
  cat(nrow(nodes.to.insert), "Nodes written\n")
  
}


DBI::dbWriteTable(con, "edgelist_intermediate_1", 
  data.frame(origin = character(0), destination = character(0),
    value = numeric(0), block_height = integer(0),
    node_index = integer(0), stringsAsFactors = FALSE), overwrite = TRUE)

base::date()
DBI::dbExecute(con, "INSERT INTO edgelist_intermediate_1 SELECT 
  origin, destination, value, block_height, node_index FROM
  edgelist JOIN nodes ON edgelist.origin = nodes.node")
base::date()


DBI::dbExecute(con, 
  "ALTER TABLE edgelist_intermediate_1 RENAME COLUMN node_index TO origin_index")


DBI::dbWriteTable(con, "edgelist_intermediate_2", 
  data.frame(origin = character(0), destination = character(0),
    origin_index = integer(0), node_index = integer(0),
    value = numeric(0), block_height = integer(0), stringsAsFactors = FALSE))

base::date()
DBI::dbExecute(con, "INSERT INTO edgelist_intermediate_2 SELECT 
  origin, destination, origin_index, node_index, value, block_height FROM
  edgelist_intermediate_1 JOIN nodes ON edgelist_intermediate_1.destination = nodes.node")
base::date()

DBI::dbExecute(con, 
  "ALTER TABLE edgelist_intermediate_2 RENAME COLUMN node_index TO destination_index")


