
catvec <- function(x) {
  # выводит на экран вектор в виде столбца (для переноса в эл. таблицу)
  xnames <- names(x)
  
  if (!is.null(xnames)) {
    for (i in 1:length(x)) {
      cat(xnames[i], "\t", x[i], "\n")
    }
  } else {
    for (z in x) {
      cat(z, "\n")
    }
  }
}
