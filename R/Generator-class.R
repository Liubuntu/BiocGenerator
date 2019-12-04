.BiocGenerator <- setRefClass(
    "BiocGenerator",
    fields = c(
        source = "ANY",
        dim = "integer",
        yieldSize = "integer",
        offset = "integer",
        totalbatch = "integer",
        nbatch = "integer"
    )
)

#' Construct, yield, and iterate over rectangular objects
#'
#' @examples
#' se <- SummarizedExperiment(matrix(rnorm(1000), 10))
#' gen <- BiocGenerator(se, yieldSize = c(NA, 20))
#' @export
BiocGenerator <-
    function(source, dim = base::dim(source), yieldSize = dim(source))
{
    yieldSize <- as.integer(yieldSize)
    offset <- rep(1L, length(dim))
    totalbatch <- ceiling(dim/yieldSize)
    totalbatch <- as.integer(prod(totalbatch[!is.na(totalbatch)]))
    .BiocGenerator(source = source, dim = dim, yieldSize = yieldSize,
                   offset = offset, totalbatch = totalbatch, nbatch = 0L)
}

.source <- function(x) x$source

.offset <- function(x) x$offset

`.offset<-` <- function(x, value) {
    updt <- !is.na(.offset(x))
    x$offset[updt] <- as.integer(value[updt])
    x
}

.totalbatch <- function(x) x$totalbatch

.nbatch <- function(x) x$nbatch

`.nbatch<-` <- function(x, value) {
    x$nbatch <- as.integer(value)
    x
}


#' @export
setMethod(
    "dim", "BiocGenerator",
    function(x)
{
    x$dim
})

## FIXME: already in Rsamtools
## setGeneric("yieldSize", function(object, ...) standardGeneric("yieldSize"))
#' @importFrom Rsamtools yieldSize
#' @export
setMethod(
    "yieldSize", "BiocGenerator",
    function(object)
{
    object$yieldSize
})

.loopseq <- function(from, to, ceiling) {
    c(seq(from, ceiling), seq_len(to))
}

## "yield" returns a chunk of SummarizedExperiment object, need to
## realize any lazy data here. 

#' @export
setGeneric("yield", function(object) standardGeneric("yield"))
setMethod(
    "yield", "BiocGenerator",
    function(object)
{
    ## browser()
    nb <- .nbatch(object) + 1L
    if (nb > .totalbatch(object))
        return(NULL)
    .nbatch(object) <- nb
    offset <- pmin(.offset(object), dim(object))
    ## FIXME: better stopping behavior
    ## if (any(offset >= dim(object), na.rm=TRUE))
    ##     return(NULL)
    ## to <- pmin(offset - 1L + yieldSize(object), dim(object))
    to <- offset - 1L + yieldSize(object)
    idx <- vector("list", length(offset))
    idx[is.na(to)] <- rep(list(TRUE), sum(is.na(to)))
    idx[!is.na(to)] <- Map(seq, offset[!is.na(to)], to[!is.na(to)])
    if (any(to >= dim(object), na.rm = TRUE)) {
        to <- to - dim(object)
        idx[!is.na(to)] <- Map(.loopseq, offset[!is.na(to)],
                               to[!is.na(to)],
                               as.list(dim(object))[!is.na(to)])
    }
    .offset(object) <- to + 1L  ## update the next offset after yielding.
    do.call("[", c(list(.source(object)), idx))
})

.pretty_dim <- function(x)
    paste(x, collapse = " ")

#' @export
setMethod(
    "show", "BiocGenerator",
    function(object)
{
    cat(
        "class: ", class(object), "\n",
        "dim: ", .pretty_dim(dim(object)), "\n",
        "yieldSize: ", .pretty_dim(yieldSize(object)), "\n",
        "current offset: ", .pretty_dim(.offset(object)), "\n",
        "total batches: ", .totalbatch(object), "\n", 
        "nbatch: ", .nbatch(object), "\n",
        sep = ""
    )
})

#' @export
setGeneric(
    "iterate",
    function(x, fun, ...) standardGeneric("iterate"),
    signature = "x"
)


#' @examples
#' iterate(gen, colSums)

setMethod(
    "iterate", "BiocGenerator",
    function(x, fun, ...)
{
    ## browser()
    result <- list()
    repeat {
        value <- yield(x)
        if (is.null(value))
            break
        result <- c(result, list(value))
    }
    result
})

##
BiocGenerator <- function()
setMethod(
    "", "BiocGenerator")

SEGenerator <- function(SE, label_column, yieldSize, genFun) {
    gen <- BiocGenerator(SE, yieldSize = yieldSize)
    function() {
        se_chunk <- yield(gen)
        x <- genFun(se_chunk)
        y <- se_chunk[[label_column]]
        if (!is(y, "numeric")) {
            y <- factor(y, levels = unique(SE[[label_column]]))
            y <- as.matrix(model.matrix(~y+0))
            colnames(y) <- sub("^y", "", colnames(y))
        }
        list(x, y)  ## return a realized X as matrix
    }
}

cpmNorm <- function(SE)
{
    dat <- counts(SE)
    t(edgeR::cpm(dat, log = TRUE))
}

## SEGEnerator(SE, "label", c(NA, 20), cpmNorm)
## ERROR in: ~/Documents/Research/BiocGenerator-all/BiocGenerator-demo.Rmd
