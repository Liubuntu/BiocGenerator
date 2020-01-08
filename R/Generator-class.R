.BiocGenerator <- setRefClass(
    "BiocGenerator",
    fields = c(
        source = "ANY",
        dim = "integer",
        yieldSize = "integer",
        offset = "integer"
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
    .BiocGenerator(source = source, dim = dim, yieldSize = yieldSize,
                   offset = offset)
}

.source <- function(x) x$source

.offset <- function(x) x$offset

`.offset<-` <- function(x, value) {
    updt <- !is.na(.offset(x))
    x$offset[updt] <- as.integer(value[updt])
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

## now "yield" is a non-stopping yielding, will loop into beginning. 

#' @export
setGeneric("yield", function(object, ...) standardGeneric("yield"))

#' @export
setMethod(
    "yield", "BiocGenerator",
    function(object, stop = FALSE)
{
    ## browser()
    offset <- .offset(object)
    if (stop && offset[!is.na(offset)] == dim(object)[!is.na(offset)])
        return(NULL)
    ## offset <- pmin(.offset(object), dim(object))
    ## FIXME: better stopping behavior
    ## if (any(offset >= dim(object), na.rm=TRUE))
    ##     return(NULL)
    ## to <- pmin(offset - 1L + yieldSize(object), dim(object))
    to <- offset - 1L + yieldSize(object)  ## to... 
    idx <- vector("list", length(offset))
    idx[is.na(to)] <- rep(list(TRUE), sum(is.na(to)))
    if (any(to >= dim(object), na.rm = TRUE)) {
        if (!stop) {
            to <- to - dim(object) ## reset "to" if >= dim
            idx[!is.na(to)] <- Map(.loopseq, offset[!is.na(to)],
                                   to[!is.na(to)],
                                   as.list(dim(object))[!is.na(to)])
            .offset(object) <- to + 1L
        } else {
            idx[!is.na(to)] <- Map(seq, offset[!is.na(to)],
                                   as.list(dim(object))[!is.na(to)])
            .offset(object)[!is.na(to)] <- dim(object)[!is.na(to)]
        }
    } else {
        idx[!is.na(to)] <- Map(seq, offset[!is.na(to)], to[!is.na(to)])
        .offset(object) <- to + 1L  ## update .offset(object) after yielding.
    }
    do.call("[", c(list(.source(object)), idx)) ## return chunk of .source(object)
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

## because "yield" is non-stopping, iterate will be non-stopping... 
setMethod(
    "iterate", "BiocGenerator",
    function(x, fun, stop = FALSE, ...) ## need to make sure "fun" is a method for
                          ## class(x).
{
    ## browser()
    result <- list()
    repeat {
        x_chunk <- yield(x, stop = stop) 
        if (is.null(x_chunk))
            break
        value <- fun(x_chunk)
        result <- c(result, list(value))
    }
})

## SEGenerator has added a "y" column, and a "genFun" based on the
## BiocGenerator ref class. It returns a function, which returns a
## list, [[1]] being the yielded (realized) data chunk by applying the
## "genFun", and [[2]] being "y" label. The function is ready to pass
## into the deep learning model.

SEGenerator <- function(SE, label_column, yieldSize, genFun) {
    ## genFun: any function, e.g., for normalization
    ## label_column: usually a categorical variable used as Y ~
    ## yieldSize: a vector indicating number of rows and column to yield each time.
    gen <- BiocGenerator(SE, yieldSize = yieldSize)
    function() {
        se_chunk <- yield(gen)  ## non-stopping, looped yielding. 
        x <- genFun(se_chunk)  ## e.g., cpm normalization
        y <- se_chunk[[label_column]]  ## categorical variable 
        if (!is(y, "numeric")) {
            y <- factor(y, levels = unique(SE[[label_column]]))
            y <- as.matrix(model.matrix(~y+0))
            colnames(y) <- sub("^y", "", colnames(y))
        }
        list(x, y)  ## return a realized X (data chunk) as matrix, and label y.
    }
}

## cpmNorm is defined directly on SummarizedExperiment object, which must include an assay name of "count"
cpmNorm <- function(SE)  
{
    dat <- counts(SE)
    t(edgeR::cpm(dat, log = TRUE))
}

## SEGEnerator(SE, "label", c(NA, 20), cpmNorm)
## ERROR in: ~/Documents/Research/BiocGenerator-all/BiocGenerator-demo.Rmd
