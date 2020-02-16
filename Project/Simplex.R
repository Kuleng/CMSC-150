
colnames = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l","m", "n", "o", "DenverTotal", "PhoenixTotal", "DallasTotal", "SacramentoTotal", "SaltLakeTotal", "Albuquerque", "Chicago", "New York")
mainMatrix = matrix(0, nrow = 1, ncol = length(colnames), dimnames = list(c(0), colnames))

print(mainMatrix)