library(data.table)
library(rbch)
library(future.apply)

is.dogecoin <- FALSE

blockchain.conf.file <- ""
# Input filepath for your {blockchain}.conf file

data.dir <- ""
# Input data directory here, with trailing "/"

current.block.height <- NA_integer_
# current.block.height <- rbch::getblockchaininfo(bch.config)@result$blocks

n.threads <- min(c(6, parallelly::availableCores()))
# Recommended no more than 6 threads since all threads query the single blockchain daemon process.


blockchain.config <- rbch::conrpc(blockchain.conf.file)
rpcport <- readLines(blockchain.conf.file)
rpcport <- rpcport[grepl("rpcport", rpcport) ]
if (length(rpcport) > 0) {
  blockchain.config@url <- paste0("http://127.0.0.1:", gsub("[^0-9]", "", rpcport))
}



getblock.doge <- function(con, blockhash, verbosity = TRUE) {
  bh <- as.character(blockhash)
  pl <- unname(list(blockhash = bh, verbosity = verbosity))
  rpcpost(con, "getblock", pl)
}


future::plan(future::multiprocess(workers = n.threads))

block.times <- future.apply::future_lapply(0:current.block.height, function(iter.block.height) {
  
  if (iter.block.height %% 1000 == 0) {
    cat(iter.block.height, base::date(), "\n")
  }
  
  block.hash <- rbch::getblockhash(blockchain.config, iter.block.height)
  
  if ( ! is.dogecoin) {
    block.data <- rbch::getblock(blockchain.config, blockhash = block.hash@result, verbosity = "l1")
  } else {
    block.data <- getblock.doge(blockchain.config, blockhash = block.hash@result, verbosity = TRUE)
  }
  data.frame(block_height = iter.block.height, block_time = block.data@result$time)
})

block.times <- data.table::rbindlist(block.times)

saveRDS(block.times, file = paste0(data.dir, "block_times.rds"))



