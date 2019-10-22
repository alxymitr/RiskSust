analyse_table_expolin <- function(da) {
    #' анализирует все показатели из таблицы (строки -- показатели, столбцы -- годы)
    #' @param da таблица данных (после фильтрации показателей из паспорта МО)
    
    dat <- t(da)
    dat_dim <- dim(dat)
    
    # анализируемая таблица
    da1 <- dat[3:dat_dim[1], ]
    
    # краткие обозначения показателей
    pokn <- paste0("v", 1:dat_dim[2])
    
    colnames(da1) <- pokn
    da1 <- cbind(year = rownames(da1), da1)
    da1 <- data.frame(apply(da1, 2, as.numeric), stringsAsFactors = FALSE)
    
    # таблица результатов
    resu <- cbind(dat[1, ], dat[2, ], NA, NA, NA, NA)
    rownames(resu) <- pokn
    colnames(resu) <- c("Показатель", "Ед. измерения", "Код выхода", 
        "Вер. роста", "Нач. прог.", "Кон. прог.")
    
    for (i in 1:length(pokn)) {
        cat("\n==========================================\n")
        cat(pokn[i], "\n")
        cat(dat[1, i], "\n")
        exle <- est_expolin(da1, pokn[i], yl = dat[2, i])
        resu[i, "Код выхода"] <- exle$e_code
        resu[i, "Вер. роста"] <- round(exle$mprob_gro, 4)
        resu[i, "Нач. прог."] <- exle$fci[1]
        resu[i, "Кон. прог."] <- exle$fci[2]
        rownames(resu)[i] <- exle$vn
    }
    resu
}
