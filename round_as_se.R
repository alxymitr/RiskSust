
round_as_se <- function(v,
                        err,
                        orc = 0,
                        ndi = 2,
                        decs_comma = TRUE) {
  # округляет действительное число
  # в соответствии со станд. ошибкой
  #
  # in:
  #   v: значение
  #   err: ст. ошибка
  #   orc: корректирующий показатель степени (будет вынесен)
  #   ndi: число значащих цифр в ст. ошибке
  #   decs_comma: десятичный разделитель ","
  
  if (orc != 0) {
    v <- v * 10^(-orc)
    err <- err * 10^(-orc)
  }
  
  err_split <- strsplit(format(err, scientific = T), "e")[[1]]
  orde <- as.integer(err_split[2])
  err1 <- signif(err, ndi)
  vr <- round(v, -orde + ndi - 1)
  ve <- paste0(vr, " ± ", err1)
  if (decs_comma) {
    ve <- gsub("\\.", ",", ve)
  }
  
  if (orc != 0) {
    ve <- paste0("(",ve,") × 10^",orc)
  }
  
  list(ve = ve, v = vr * 10^orc, err = err1 * 10^orc)
}
