
library(data.table)
library(RSQLite)
library(DBI)
# NOTE: Also need lubridate package installed, but not loading it due to 
# it masking functions

is.btc <- FALSE
# Change to TRUE if processing BTC

data.dir <- ""
# Input data directory here, with trailing "/"


con <- DBI::dbConnect(RSQLite::SQLite(), paste0(data.dir, "tx-graph-node-indices.db"))


if ( ! is.btc) {
  
  master.edgelist.output.created <- master.edgelist[
    (! is.na(master.edgelist$value)) & master.edgelist$value > 0 , c("destination_index", "block_height")]
  colnames(master.edgelist.output.created) <- c("output_index", "output.created.block_height")
  setDT(master.edgelist.output.created)
  
  master.edgelist.output.spent <- master.edgelist[, c("origin_index", "block_height")]
  colnames(master.edgelist.output.spent) <- c("output_index", "output.spent.block_height")
  setDT(master.edgelist.output.spent)
  
  # Only include positive _value_s for output created, since that's the value of the created output
  # Then, below do an "inner merge" to get the proper outputs on the spending side
  
  rm(master.edgelist)
  master.edgelist.output.spent <- merge(master.edgelist.output.created, master.edgelist.output.spent)
  rm(master.edgelist.output.created)
  
} else {
  
  # Split data requests and processing if BTC to avoid integer overflow issue in this step of processing.
  # data.table issue:
  # https://github.com/Rdatatable/data.table/issues/3957
  # SQLite issue: "negative length vectors are not allowed"
  
  master.edgelist.lt.500k <- DBI::dbGetQuery(con, 
    "SELECT origin_index, destination_index,block_height,value FROM edgelist_intermediate_2 WHERE block_height < 500000")
  
  master.edgelist.output.created.lt.500k <- master.edgelist.lt.500k[
    (! is.na(master.edgelist.lt.500k$value)) & master.edgelist.lt.500k$value > 0 , c("destination_index", "block_height")]
  colnames(master.edgelist.output.created.lt.500k) <- c("output_index", "output.created.block_height")
  setDT(master.edgelist.output.created.lt.500k)
  
  master.edgelist.output.spent.lt.500k <- master.edgelist.lt.500k[, c("origin_index", "block_height")]
  colnames(master.edgelist.output.spent.lt.500k) <- c("output_index", "output.spent.block_height")
  setDT(master.edgelist.output.spent.lt.500k)
  
  rm(master.edgelist.lt.500k)
  
  
  master.edgelist.gt.500k.lt.700k <- DBI::dbGetQuery(con, 
    "SELECT origin_index, destination_index,block_height,value FROM edgelist_intermediate_2 WHERE block_height > 499999 AND block_height < 700000")
  
  master.edgelist.output.created.gt.500k.lt.700k <- master.edgelist.gt.500k.lt.700k[
    (! is.na(master.edgelist.gt.500k.lt.700k$value) & master.edgelist.gt.500k.lt.700k$value > 0 ) , c("destination_index", "block_height")]
  colnames(master.edgelist.output.created.gt.500k.lt.700k) <- c("output_index", "output.created.block_height")
  setDT(master.edgelist.output.created.gt.500k.lt.700k)
  
  master.edgelist.output.spent.gt.500k.lt.700k <- master.edgelist.gt.500k.lt.700k[, c("origin_index", "block_height")]
  colnames(master.edgelist.output.spent.gt.500k.lt.700k) <- c("output_index", "output.spent.block_height")
  setDT(master.edgelist.output.spent.gt.500k.lt.700k)
  
  rm(master.edgelist.gt.500k.lt.700k)
  
  
  master.edgelist.gt.700k <- DBI::dbGetQuery(con, 
    "SELECT origin_index, destination_index,block_height,value FROM edgelist_intermediate_2 WHERE block_height > 699999")
  
  master.edgelist.output.created.gt.700k <- master.edgelist.gt.700k[
    (! is.na(master.edgelist.gt.700k$value) & master.edgelist.gt.700k$value > 0 ), c("destination_index", "block_height")]
  colnames(master.edgelist.output.created.gt.700k) <- c("output_index", "output.created.block_height")
  setDT(master.edgelist.output.created.gt.700k)
  
  master.edgelist.output.spent.gt.700k <- master.edgelist.gt.700k[, c("origin_index", "block_height")]
  colnames(master.edgelist.output.spent.gt.700k) <- c("output_index", "output.spent.block_height")
  setDT(master.edgelist.output.spent.gt.700k)
  
  rm(master.edgelist.gt.700k)
  
  
  master.edgelist.output.created <- rbindlist(list(master.edgelist.output.created.lt.500k, 
    master.edgelist.output.created.gt.500k.lt.700k, master.edgelist.output.created.gt.700k))
  
  master.edgelist.output.spent.lt.500k <- merge(master.edgelist.output.created, master.edgelist.output.spent.lt.500k)
  master.edgelist.output.spent.gt.500k.lt.700k <- merge(master.edgelist.output.created, master.edgelist.output.spent.gt.500k.lt.700k)
  master.edgelist.output.spent.gt.700k <- merge(master.edgelist.output.created, master.edgelist.output.spent.gt.700k)
  
  master.edgelist.output.spent <- rbindlist(list(master.edgelist.output.spent.lt.500k, 
    master.edgelist.output.spent.gt.500k.lt.700k, master.edgelist.output.spent.gt.700k))
  
  rm(master.edgelist.output.created)
  
}



block.times <- readRDS(paste0(data.dir, "block_times.rds"))
colnames(block.times) <- c("output.created.block_height", "output.created.block_time")
setDT(block.times)

master.edgelist.output.spent <- merge(master.edgelist.output.spent, block.times, by = "output.created.block_height")

colnames(block.times) <- c("output.spent.block_height", "output.spent.block_time")

master.edgelist.output.spent <- merge(master.edgelist.output.spent, block.times, by = "output.spent.block_height")

# TODO: Explore phenomenon of out-of-order block timestamps and decide what to do about them

master.edgelist.output.spent[, output.spend.age := output.spent.block_time - output.created.block_time]

# These reduce RAM usage if desired:
# master.edgelist.output.spent[, output.created.block_time := NULL]
# master.edgelist.output.spent[, output.created.block_height := NULL]
# master.edgelist.output.spent[, output_index := NULL]


output.spent.block_time.intermediate <- unique(master.edgelist.output.spent[, .(output.spent.block_time)])

output.spent.block_time.intermediate[, output.spent.block_time.week :=
    lubridate::isoweek(as.POSIXct(output.spent.block_time, origin = "1970-01-01", tz = "UTC"))]

output.spent.block_time.intermediate[, output.spent.block_time.isoweekyear :=
    lubridate::isoyear(as.POSIXct(output.spent.block_time, origin = "1970-01-01", tz = "UTC"))]

output.spent.block_time.intermediate[,
  output.spent.block_time.week := factor(paste0(output.spent.block_time.isoweekyear, "-", 
    formatC(output.spent.block_time.week, width = 2, flag = "0")))]

master.edgelist.output.spent <- merge(master.edgelist.output.spent, 
  output.spent.block_time.intermediate[, .(output.spent.block_time, output.spent.block_time.week)], by = "output.spent.block_time")

saveRDS(master.edgelist.output.spent, paste0(data.dir, "master_edgelist_output_spent.rds"))

saveRDS(master.edgelist.output.spent, paste0(data.dir, "master_edgelist_output_spent-uncompressed.rds"), compress = FALSE)


