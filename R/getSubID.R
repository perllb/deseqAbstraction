#' @name getSubID
#' @description takes a vector of RepeatMasker TE_IDs with unique _x tag, and removes the tag.
#' @param uniqID: vector with unique IDs
#' @title getSubID - easy parsing
#' @export getSubID
#' @examples
#' getSubID(uniqID = rownames(countTEdata))
 
# get the unique name f
getSubID = function(uniqID=NULL){
  return(sub("_[^_]+$", "", uniqID))
}