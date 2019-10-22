
#' оценивает параметры (4) эксполинейной функции, прогнозирует, находит производную
#' для будущих значений для случайных оценок параметров
#' оценивает вероятность роста прогноза 
#' разработчик А.Ю. Митрофанов
#'
#' История модификаций
#' -------------------
#' v.0.6 в случае, когда анализ не удался (ранее e_code=-2),
#'       сглаживает ряд фильтром Ходрика-Прескотта и запускается снова
#'       e_code=-5: сглаживание не удалось (есть NA в середине ряда)
#'		 e_code=-4: не удалось проанализировать сглаженный ряд
#'       выдаёт название переменной ($vn)
#'       результат выдаётся видимо
#'       выдаёт прогнозный интервал ($fci)
#'
#' v.0.5 возвращается список (e_code, mprob_gro)
#'       e_code -- код выхода, mprob_gro -- средняя вероятность роста
#'       e_code= 0: наилучший результат
#'		 e_code=-1: возможно, эксполинейная модель не подходит (широкие границы для производной)
#'		 e_code=-2: NA в значениях производной для случайных наборов параметров
#'                  эксполин. тренд не подходит
#'       e_code=-3: в значениях показателя нет вариации, mprob_gro = 0.5
#'
#' @param df таблица данных со столбцами year, vy
#' @param vy зависимая переменная (строка)
#' @param ycast горизонт прогнозирования (лет)
#' @param yl этикетка оси ординат на графике
#' @param legpos положение легенды на графике
#' @param nrand число случайных стартовых значений (оценка параметров)
#' @param nb_rand_parameters число случайных наборов параметров (прогнозные интервалы, производная)
#' @param nb_year_subi число подынтервалов, на которые делится каждый прогнозный год (вычисление производной)
#' @param fc_omega надёжность прогнозных интервалов
#' @param ra_level граница (%), при превышении которой вверх или вниз(-) выдаётся предупреждение о ненадёжности модели
#' @param orig_series индикатор того, что анализируется оригинальный ряд
#'                    (не сглаженный фильтром Ходрика-Прескотта)
#'
#' @return
#' @export
#'
#' @examples
est_expolin <- function(df, vy,
	                    ycast = 5,
	                    yl = vy,
	                    legpos = "bottomright",
	                    nrand = 1000,
                        nb_rand_parameters = 1000,
                        nb_year_subi = 4,
                        fc_omega = 0.9,
                        ra_level = 50,
                        orig_series = TRUE) {
    
    expol <- function(b, x) {
        # значения эксполинейной функции на векторе x
        b[1] + b[2] * (x - b[4])/(1 + exp(-b[3] * (x - b[4])))
    }
    
    expol_deriv <- function(b, x) {
        # значения производной эксполинейной функции на векторе x
        z <- exp(b[3] * (b[4] - x))
        b[2]/(z + 1) * (1 - b[3] * (b[4] - x) * z/(z + 1))
    }
    
    rss <- function(b, x, y) {
        # сумма квадратов остатков
        yhat <- expol(b, x)
        1000 * sum((y - yhat)^2)/length(yhat)
    }
    
    gradi <- function(x, b) {
        # градиент по b в точке x
        z <- exp(b[3] * (b[4] - x))
        z1 <- z + 1
        z12 <- z/z1^2
        b4x <- b[4] - x
        c(1, -b4x/z1, b[2] * b4x^2 * z12, b[2] * (b[3] * b4x * z12 - 1/z1))
    }
    
    #' таблица оценок коэффициентов, ст. ошибок
    #'
    #' @param vv имена предикторов
    #' @param b оценки коэфф-тов
    #' @param seb ст. ошибки коэфф-тов
    #' @param dfR ЧСС остатков
    #'
    #' @return
    #' @export
    #'
    #' @examples
    summary_table <- function(vv, b, seb, dfR, rus = FALSE) {
        df <- data.frame(b = b, seb = seb, tstat = b/seb, pv = 2 * (1 - pt(abs(b/seb), 
            df = dfR)))
        
        rownames(df) <- vv
        if (rus) {
            colnames(df) <- c("Оценка", "Ст.Ошибка", "t", "Pr(>|t|)")
            cat("Коэффициент ± ст. ошибка:\n")
        } else {
            colnames(df) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
            cat("Coefficient ± st. error:\n")
        }
        
        for (i in 1:length(b)) {
            coese <- round_as_se(df[i, 1], df[i, 2])$ve
            cat(vv[i], ":", coese, "\n")
        }
        cat("\n")
        print(signif(df, 5))
    }
    
	hp_smooth <- function(df, vn) {
	    #' сглаживает ряд (годовой) фильтром Ходрика-Прескотта
	    #'
	    #' @param df таблица данных со столбцами year, vn
	    #' @param vn зависимая переменная 
	    #'    
	    #' @return если в vn нет 'дыр' (NAs допускаются только в начале и в конце),
	    #'         таблица данных со столбцами year, vn, vn с добавкой '_hp'
	    #'         иначе NULL
	    
	    df_names <- colnames(df)
	    
	    stopifnot("year" %in% df_names)
	    stopifnot(vn %in% df_names)

	    df0 <- df[, c("year", vn)]
	    
	    df1 <- df0[complete.cases(df0), ]
	    n <- nrow(df1)
	    
	    if (!all.equal(df1$year, df1$year[1]:df1$year[n])) {
	        # есть 'дыры' -- NAs
	        cat("\n*** Есть отсутствующие значения (NA) в середине ряда. ***\n")
	        cat("*** сглаживание отменено ***\n")
	        return(NULL)
	    }
	    
	    vn_ts <- ts(df1[, vn], start = df1$year[1], freq = 1)
	    vn_ts_hp <- mFilter::hpfilter(vn_ts, freq = 1)
	    # plot(vn_ts_hp)
	    
	    smo <- as.vector(vn_ts_hp$trend)
	    
	    plot1(df1[, vn] ~ df1[, "year"], type = "o", lwd = 2, pch = 19, xlab = "Год", ylab = vn)
	    lines1(smo ~ df1[, "year"], type = "o", lty = 2)
	    grid()
	    legend1("bottomright", legend = c("Наблюдаемый ряд", "Сглаженный (H-P) ряд"), 
	        lty = c(1, 2), pch = c(19, 1), lwd = 2)
	    
	    df0[complete.cases(df0), paste0(vn, "_hp")] <- smo
	    df0
	}

    cat("+-----------------------------------------+\n")
    cat("| Версия 0.6, разработчик А.Ю. Митрофанов |\n")
    cat("+-----------------------------------------+\n")
    
    cat("Анализируемый ряд:", vy, "\n\n")

    # проверка наличия столбцов 'year', vy
    dfvars <- colnames(df)
    stopifnot("year" %in% dfvars)
    stopifnot(vy %in% dfvars)
    
    df <- df[, c("year", vy)]

    # сохранение для возможного последующего сглаживания
    df0 <- df

    df <- df[complete.cases(df), ]
    colnames(df)[2] <- "y"

    # прогнозный интервал
    fci <- max(df$year) + c(1, ycast)
    
    if (length(unique(df$y)) == 1) {
    	cat("*** В значениях показателя нет вариации ***\n")
    	return(list(e_code = -3, mprob_gro = 0.5, vn = vy, fci = fci))
    }

    # нормировка значений
    yeara <- min(df$year)
    yearb <- max(df$year)
    ya <- min(df$y)
    yb <- max(df$y)
    
    x <- (df$year - yeara)/(yearb - yeara)
    y <- (df$y - ya)/(yb - ya)
    
    # оценивание параметров (в нормированной шкале) для случайных стартовых значений:
    # 4 стартовых, 4 оптимизированных, сходимость, оптимальное значение
    resu <- matrix(NA, nrand, 10)
    
    set.seed(1234)
    resu[, 1:4] <- matrix(runif(nrand * 4, min = -5, max = 5), ncol = 4)
    
    for (i in 1:nrand) {
        b <- resu[i, 1:4]
        opt <- optim(b, rss, gr = NULL, x, y, control = list(maxit = 1000))
        resu[i, ] <- c(b, opt$par, opt$convergence, opt$value)
    }
    
    # выделяются только те варианты стартовых значений, в которых достигнута
    # сходимость
    resu <- resu[!resu[, 9], ]
    stopifnot(nrow(resu) > 0)
    
    # случайный набор параметров, для которого RSS минимальна
    b <- resu[which.min(resu[, 10]), 5:8]
    
    cat("Использовано", nrand, "случайных векторов стартовых значений.\n")
    
    # матрица градиентов эксполин. функции в точках наблюдений (аналог матрицы плана)
    xlen <- length(x)
    X <- matrix(NA, nrow = xlen, ncol = 4)
    for (i in 1:xlen) {
        X[i, ] <- gradi(x[i], b)
    }
    
    # оценка дисперсии случ. возмущения
    dfR <- xlen - 4
    RSS <- rss(b, x, y)/1000 * xlen
    sigma2hat <- RSS/dfR
    sigmahat <- sqrt(sigma2hat)
    
    # ковариационная матрица и ст. ошибки оценок коэффициентов
	varb <- sigma2hat * chol2inv(chol(t(X) %*% X))

    seb <- sqrt(diag(varb))
    TSS <- sum((y - mean(y))^2)
    
    if (FALSE) {
        cat("Параметры в нормированной шкале\n")
        # print(summary(lm(y ~ x)))
        summary_table(paste0("b", 0:3), b, seb, dfR, rus = TRUE)
        cat("dfR=    ", dfR, "\n")
        cat("sigma^= ", signif(sigmahat, 5), "\n")
        cat("R^2=    ", round(1 - RSS/TSS, 4), "\n")
        cat("R^2adj.=", round(1 - sigma2hat/var(y), 4), "\n")
        cat("\n*** Все величины приблизительны, поскольку модель\n")
        cat("не учитывает автокоррелированность случайных возмущений. ***\n")
    }
    
    #----- случайные параметры (четвёрки в строках) в нормированной шкале ----------
    
    b_rand <- mvtnorm::rmvnorm(nb_rand_parameters, mean = b, sigma = varb)
    cat("Использовано", nb_rand_parameters, "случайных векторов параметров.\n")
    
    # прогнозные значения времени в норм. шкале (мелкий шаг, для производной)
    xdcast <- seq(from = max(df$year), to = max(df$year) + ycast, by = 1/nb_year_subi)
    xdcast1 <- (xdcast - yeara)/(yearb - yeara)
    
    # график в нормированной шкале
    if (FALSE) {
        plot(y ~ x, type = "n", ylim = c(-0.5, 1.5), xlim = c(x[1], xdcast1[length(xdcast1)]))
        yhat <- expol(b, xdcast1)
        
        for (i in 1:nb_rand_parameters) {
            yhat <- expol(b_rand[i, ], xdcast1)
            lines(yhat ~ xdcast1, type = "l", col = "grey")
        }
        
        lines(y ~ x, type = "o")
        yhat <- expol(b, x)
        lines(yhat ~ x, type = "o", col = "red")
    }
    
    #----- прогнозные значения производной нормированной эксполинейной функции --------
    
    # при случайных значениях параметров
    der_rand <- matrix(NA, nb_rand_parameters, length(xdcast))
    for (i in 1:nb_rand_parameters) {
        der_rand[i, ] <- expol_deriv(b_rand[i, ], xdcast1)
    }
    
    if (any(is.nan(der_rand))) {
        cat("\n-------------------------------------------------------------\n")
        cat("В значениях производной при случайных значениях параметров\n")
        cat("имеются отсутствующие (NA) => эксполинейный тренд неприменим.\n")
        cat("-------------------------------------------------------------\n")

        if (orig_series) {
        	# ряд оригинальный
        	cat("\n*** Сглаживание ряда фильтром Ходрика-Прескотта ***\n\n")
        	df_hp <- hp_smooth(df0, vy)
        	if (!is.null(df_hp)) {
        		vy_new <- paste0(vy, "_hp")	
        		resu_hp <- est_expolin(df_hp, vy_new, ycast, yl = vy_new,
        			legpos, nrand, nb_rand_parameters, nb_year_subi,
        			fc_omega, ra_level, orig_series = FALSE)
        		return(resu_hp)
        	} else {
        		# сглаживание не удалось (есть NA в середине ряда)
        		return(list(e_code = -5, mprob_gro = NA, vn = vy, fci = fci))
        	}
        } else {
        	# ряд сглаженный -- выходить
        	return(list(e_code = -4, mprob_gro = NA, vn = vy, fci = fci))
        }
    }
    
    # статистика производных на частой (нормированной) шкале
    der_rand_med <- apply(der_rand, 2, median)
    
    der_rand_q25 <- apply(der_rand, 2, function(x) {
        quantile(x, probs = 0.25)
    })
    
    der_rand_q75 <- apply(der_rand, 2, function(x) {
        quantile(x, probs = 0.75)
    })
    
    par(mfrow = c(3, 1))
    
    # график медианных значений производной с квартилями (в нормированной шкале)
    plot1(der_rand_med ~ xdcast, type = "l", ylim = c(min(der_rand_q25), max(der_rand_q75)), 
        xlab = "Год", ylab = "Производная", lwd = 2, main = paste0(vy, 
            ": ", "прогноз производной (в норм. шкале)"))
    lines1(der_rand_q25 ~ xdcast, lty = 2)
    lines1(der_rand_q75 ~ xdcast, lty = 2)
    grid()
    abline(h = 0)
    legend1("topright", legend = c("Медиана", "Квартили"), lty = c(1, 
        2), lwd = 2)
    
    # индикатор неотрицательной производной
    der_rand01 <- (der_rand >= 0) * 1
    
    # вероятность роста
    der_rand01up <- apply(der_rand01, 2, mean)
    der_rand01up_mean <- mean(der_rand01up)
    
    # график динамики вероятности роста
    plot1(der_rand01up ~ xdcast, type = "l", ylim = c(0, 1), xlab = "Год", ylab = "Вероятность", 
        lwd = 2, main = paste0(vy, ": ", "прогноз вероятности роста"))
    grid()
    
    #----- параметры в исходной шкале ----------------------------------------
    
    cat("\nКоэффициент:\n")
    b0 <- ya + (yb - ya) * b[1]
    b1 <- (yb - ya)/(yearb - yeara) * b[2]
    b2 <- b[3]/(yearb - yeara)
    b3 <- yeara + (yearb - yeara) * b[4]
    cat("b0 :", b0, "\n")
    cat("b1 :", b1, "\n")
    cat("b2 :", b2, "\n")
    cat("b3 :", b3, "\n\n")
    b_init <- c(b0, b1, b2, b3)
    
    # ковар. матрица в исходной шкале
    dia <- c(yb - ya, (yb - ya)/(yearb - yeara), 1/(yearb - yeara), (yearb - yeara))
    varb_init <- diag(dia) %*% varb %*% diag(dia)
    # print(varb_init)
    
    # ст. ошибки в исходной шкале
    seb_init <- dia * seb
    
    summary_table(paste0("b", 0:3), b_init, seb_init, dfR, rus = TRUE)
    cat("dfR=    ", dfR, "\n")
    cat("sigma^= ", signif(sigmahat * (yb - ya), 5), "\n")
    cat("R^2=    ", round(1 - RSS/TSS, 4), "\n")
    cat("R^2adj.=", round(1 - sigma2hat/var(y), 4), "\n")
    cat("\n*** Все величины приблизительны, поскольку модель\n")
    cat("не учитывает автокоррелированность случайных возмущений. ***\n")
    
    # проверка результатов
    if (FALSE) {
        m1 <- nls(y ~ b0 + (b1 * (year - b3))/(1 + exp(-b2 * (year - b3))), data = df, 
            start = list(b0 = b0, b1 = b1, b2 = b2, b3 = b3), algorithm = "port", 
            nls.control(maxiter = 1000, warnOnly = TRUE))
        
        print(summary(m1))
        # print(coef(m1))
    }
    
    #----- график в исходной шкале ------------------------------------------
    
    # нормированные значения x, включая будущие (годовые данные)
    xcast <- c(df$year, max(df$year) + (1:ycast))
    x1cast <- (xcast - yeara)/(yearb - yeara)
    
    y1cast <- expol(b, x1cast)
    ycast <- ya + (yb - ya) * y1cast
    
    # прогнозные интервалы при случайных значениях параметров
    exl_rand <- matrix(NA, nb_rand_parameters, length(x1cast))
    for (i in 1:nb_rand_parameters) {
        exl_rand[i, ] <- expol(b_rand[i, ], x1cast)
    }
    
    exl_rand_lo <- apply(exl_rand, 2, function(x) {
        quantile(x, probs = (1 - fc_omega)/2)
    })
    
    exl_rand_hi <- apply(exl_rand, 2, function(x) {
        quantile(x, probs = (1 + fc_omega)/2)
    })
    
    exl_rand_lo <- ya + (yb - ya) * exl_rand_lo
    exl_rand_hi <- ya + (yb - ya) * exl_rand_hi
    
    lo_ra <- range((exl_rand_lo/ycast - 1) * 100)
    hi_ra <- range((exl_rand_hi/ycast - 1) * 100)
    
    cat("\nДиапазон отклонений верхней границы тренда:", 
        signif(hi_ra, 3), "%\n")
    cat("Диапазон отклонений нижней границы тренда: ", 
        signif(lo_ra, 3), "%\n")
    
    exit_code <- 0
    
    if (max(hi_ra) > ra_level) {
        cat("---------------------------------------------------------------\n")
        cat("Верхняя граница тренда отклоняется вверх более, чем на", 
            ra_level, "% !!!\n")
        cat("Модель эксполинейного тренда может быть ненадёжной!\n")
        cat("---------------------------------------------------------------\n")
        exit_code <- -1
    }
    
    if (min(lo_ra) < -ra_level) {
        cat("-------------------------------------------------------------\n")
        cat("Нижняя граница тренда отклоняется вниз более, чем на", 
            ra_level, "% !!!\n")
        cat("Модель эксполинейного тренда может быть ненадёжной!\n")
        cat("-------------------------------------------------------------\n")
        exit_code <- -1
    }
    
    cat("\n-------------------------------------\n")
    cat(vy, ": средняя прогнозная вероятность:\n")
    cat("- роста:   ", der_rand01up_mean, "\n")
    cat("- снижения:", 1 - der_rand01up_mean, "\n")
    cat("-------------------------------------\n")
    
    plot1(c(df$year, xcast), c(df$y, ycast), type = "n", xlab = "Год", ylab = yl, 
        main = paste0(vy, ": ", "наблюдения и прогноз значений показателя"))
    lines1(y ~ year, type = "o", pch = 19, data = df)
    lines1(ycast ~ xcast, type = "o", lty = 2)
    lines1(exl_rand_lo ~ xcast, lty = 3)
    lines1(exl_rand_hi ~ xcast, lty = 3)
    
    grid()
    legend1(legpos, legend = c("Наблюдаемый ряд", "Эксполинейный тренд", 
        paste0("Границы для тренда (", fc_omega, ")")), lty = c(1, 
        2, 3), pch = c(19, 1, NA), lwd = 2)
    
    return(list(e_code = exit_code, mprob_gro = der_rand01up_mean, vn = vy, fci = fci))
}
