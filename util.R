

#' convertAdjMatrixToSortedRankTSV
#' is formatted as [GENE_X][tab][GENE_Y][tab][RANK]
#' @param inputFilename 
#' @param inputFile name or full path to the file containing the matrix with the adjancency scores
#' @param outputFilename name or full path to the outputfile
#' @param desc 
#'
#' @return
#' @export
#'
#' @examples
convertAdjMatrixToSortedRankTSV <- function(inputFilename=NULL,
                                            inputFile=NULL,
                                            outputFilename=NULL,
                                            desc=TRUE){
	if(is.null(inputFile)){
		# Read the table from file
		tbl <- read.table(inputFilename)
	}else{
		tbl <- inputFile
	}
	# First column -> repeat the predictors
	firstCol  <- rep(rownames(tbl),length(colnames(tbl)))
	# Second column -> repeat the targets
	secondCol <- c()
	for ( i in 1:length(colnames(tbl))) {
		secondCol <- append(secondCol,rep(colnames(tbl)[i],length(rownames(tbl))))
	}
	thirdCol <- as.matrix(tbl)
	thirdCol <- as.vector(thirdCol)
	# Gets the indices from a desc sort on the adjancy measures
	if (desc){
		or <- sort.list(-thirdCol)
	}
	else{
		or <- sort.list(thirdCol)
	}
	# Glue everything together
	result <- cbind(firstCol,secondCol,thirdCol)
	# Convert to dataframe
	result <- as.data.frame(result)
	# Sort it using the indices obtained before
	result <-  result[or,,drop=FALSE]
	# Write to file if filename is given
	if (!is.null(outputFilename)){
		write.table(result,file=outputFilename,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
	}
	# else write to function output
	else{ 
		return(result)
	}
}

#' convertSortedRankTSVToAdjMatrix
#'
#' @param inputFilename name or full path to the TSV file
#' @param input 
#' @param outputFilename name or full path to the outputfile
#'
#' @return
#' @export
#'
#' @examples
convertSortedRankTSVToAdjMatrix <- function(inputFilename=NULL,
                                            input=NULL,
                                            outputFilename=NULL){
	# Read the table from file
	if(!is.null(inputFilename)){
		tbl <- read.table(inputFilename)
	}
	else{
		tbl <- input
	}
	# Sort by second column
	tbl <- tbl[sort.list(tbl[,2]),]	
	# Get the number of unique elements in the predictors
	pred <- unique(tbl[,1])
	# Get the number of unique elements in the targets
	targ <- unique(tbl[,2])
	# Pre allocate return matrix
	m <- matrix(0.0, length(pred),length(targ))
	# Set the col/row-names
	colnames(m) <- targ
	rownames(m) <- pred
	# Get the duplicates
	dups <- duplicated(tbl[,2])
	# Get the starIndices of another column
	startIndices <- which(FALSE == dups)
	for (i in 1:(length(startIndices)-1)){
		targetToAdd <- tbl[startIndices[i],2]
		predToAdd   <- tbl[startIndices[i]:(startIndices[i+1]-1),1]
		valuesToAdd <- tbl[startIndices[i]:(startIndices[i+1]-1),3]
		colIndex   <- which(colnames(m) %in% targetToAdd)
		tmp <- c()
		for (i in predToAdd){
			tmp <- c(tmp, which(rownames(m) == i))
		}
		rowIndexes <- tmp
		m[rowIndexes,colIndex] <- valuesToAdd
	}
	targetToAdd <- tbl[startIndices[length(startIndices)],2]
	predToAdd   <- tbl[startIndices[length(startIndices)]:length(tbl[,1]),1]
	valuesToAdd <- tbl[startIndices[length(startIndices)]:length(tbl[,1]),3]
	tmp <- c()
	for (i in predToAdd){
		tmp <- c(tmp, which(rownames(m) == i))
	}
	rowIndexes <- tmp
	colIndex   <- which(colnames(m) %in% targetToAdd)	
	m[rowIndexes,colIndex] <- valuesToAdd	
	if(!is.null(outputFilename)){
		write.table(m,outputFilename)
	}
	else{
		return(m)
	}
}

#' isPercentage
#'
#' @param number 
#'
#' @return
#' @export
#'
#' @examples
isPercentage <- function(number){
	if (!is.numeric(number) || number <=0 || number > 100 || length(number) != 1){
		return (FALSE)
	}
	else{
		return (TRUE)
	}
}

# Utility functions -------------------------------------------------------

#' isWholeNumber Check for integer number
#'
#' @param x 
#' @param tol 
#'
#' @return
#' @export
#'
#' @examples
isWholeNumber <- function(x, tol = .Machine$double.eps^0.5){
	return (is.numeric(x) &&abs(x - round(x)) < tol)
}

#' isIntegerVector Checks if a vector is entirely composed of integers
#'
#' @param integerVector 
#'
#' @return
#' @export
#'
#' @examples
isIntegerVector <- function(integerVector){
	return (!any(!isWholeNumber(integerVector)))
}

#' isPositiveIntegerVector Checks if a vector is entirely composed of positive integers
#'
#' @param integerVector 
#'
#' @return
#' @export
#'
#' @examples
isPositiveIntegerVector <- function(integerVector){
	return (isIntegerVector(integerVector) && !any(integerVector <0))
}

#' isPositiveNonZeroIntegerVector Checks if a vector is entirely composed of positive non zero integers
#'
#' @param integerVector 
#'
#' @return
#' @export
#'
#' @examples
isPositiveNonZeroIntegerVector <- function(integerVector){
	return (isIntegerVector(integerVector) && !any(integerVector <=0))
}

#' sortMatrix Sorts a matrix according to col and rownames
#'
#' @param mat 
#'
#' @return
#' @export
#'
#' @examples
sortMatrix <- function(mat){
	mat <- mat[sort.list(rownames(mat)),sort.list(colnames(mat))  ]
	return(mat)
}

#' isDirectory
#'
#' @param directory 
#'
#' @return
#' @export
#'
#' @examples
isDirectory <- function(directory){
	return(file.info(directory)$isdir)
}

#' areBothNullOrNonNull
#'
#' @param a 
#' @param b 
#'
#' @return
#' @export
#'
#' @examples
areBothNullOrNonNull <- function(a,b){
	return (is.null(a) && is.null(b) || !is.null(a) && !is.null(b) )
}

#' areAllNullOrNonNull
#'
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
areAllNullOrNonNull <- function(...){
	e <- ( (length((which(c(...)==FALSE)))))
	return (e!=1)
}

#' factorToNumber
#'
#' @param vec 
#'
#' @return
#' @export
#'
#' @examples
factorToNumber <- function(vec){
	return(as.numeric(levels(vec))[vec])
}

#' numericMatrix
#'
#' @param a 
#'
#' @return
#' @export
#'
#' @examples
numericMatrix <- function(a){
	a<- apply(a,c(1,2),function(x){as.numeric(x)})
}
