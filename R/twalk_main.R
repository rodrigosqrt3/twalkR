#' Executa o algoritmo MCMC t-walk
#'
#' Esta funcao implementa o algoritmo t-walk de Christen & Fox (2010),
#' um amostrador MCMC de proposito geral que nao requer ajuste manual.
#' A funcao pode executar multiplas cadeias MCMC independentes em paralelo
#' para acelerar a execucao e facilitar diagnosticos de convergencia.
#'
#' @param log_posterior Uma funcao que recebe um vetor de parametros como
#'   primeiro argumento e retorna o log da densidade posterior (escalar).
#'   Argumentos adicionais podem ser passados para esta funcao via `...`.
#' @param n_iter O numero de iteracoes a serem executadas para cada cadeia.
#' @param x0 Um vetor numerico com os valores iniciais para o primeiro ponto (`x`).
#' @param xp0 Um vetor numerico com os valores iniciais para o segundo ponto (`x'`).
#' @param n_cadeias O numero de cadeias MCMC independentes a serem executadas.
#'   O padrao e `1`, que executa uma unica cadeia em modo sequencial. Se for
#'   maior que 1, ativa o modo paralelo.
#' @param n_cores O numero de nucleos de CPU a serem usados no modo paralelo.
#'   Se `NULL` (padrao), tentara usar todos os nucleos disponiveis menos um.
#' @param ... Argumentos adicionais a serem passados para a funcao `log_posterior`.
#'
#' @return Uma lista contendo:
#' \item{todas_amostras}{Uma matriz com as amostras combinadas de todas as cadeias.}
#' \item{taxa_aceitacao}{A taxa de aceitacao media entre todas as cadeias.}
#' \item{num_iteracoes_total}{O numero total de amostras geradas (n_iter * n_cadeias).}
#' \item{dimensao}{A dimensao do espaco de parametros.}
#' \item{cadeias_individuais}{Se `n_cadeias > 1`, uma lista contendo os resultados
#'       brutos de cada cadeia separadamente, util para diagnosticos como o R-hat.}
#'
#' @export
#' @importFrom parallel detectCores makeCluster clusterEvalQ clusterExport parLapply stopCluster
#' @importFrom stats rnorm runif
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @examples
#' # Exemplo 1: Amostrar de uma Normal Bivariada (modo sequencial)
#' # E necessario o pacote 'mvtnorm' para este exemplo
#' if (requireNamespace("mvtnorm", quietly = TRUE)) {
#'   log_post <- function(x) {
#'     mvtnorm::dmvnorm(x, mean = c(0, 0), sigma = matrix(c(1, 0.8, 0.8, 1), 2, 2), log = TRUE)
#'   }
#'
#'   # Rodar com menos iteracoes para um exemplo rapido
#'   resultado_seq <- twalk(log_posterior = log_post, n_iter = 5000,
#'                          x0 = c(-1, 1), xp0 = c(1, -1))
#'
#'   plot(resultado_seq$todas_amostras, pch = '.', main = "Amostras t-walk (Sequencial)")
#' }
#'
#' \dontrun{
#' # Exemplo 2: O mesmo problema em paralelo (rodara mais rapido)
#' # Usando 2 cadeias. n_iter agora e por cadeia.
#' if (requireNamespace("mvtnorm", quietly = TRUE)) {
#'   resultado_par <- twalk(log_posterior = log_post, n_iter = 2500,
#'                          x0 = c(-1, 1), xp0 = c(1, -1), n_cadeias = 2)
#'
#'   plot(resultado_par$todas_amostras, pch = '.', main = "Amostras t-walk (Paralelo)")
#' }
#' }
twalk <- function(log_posterior, n_iter, x0, xp0,
                  n_cadeias = 1, n_cores = NULL, ...) {

  # Capturar todos os argumentos extras em uma lista
  argumentos_extras <- list(...)

  # --- BLOCO SEQUENCIAL ---
  if (n_cadeias == 1) {

    is_internal_call <- "internal_call" %in% names(argumentos_extras)

    if (!is_internal_call) {
      cat("--- Executando T-Walk em modo sequencial (1 cadeia) ---\n")
    }

    dimensao <- length(x0)

    # Cria uma cópia limpa dos argumentos extras para uso interno,
    # removendo o sinalizador 'internal_call'.
    argumentos_internos <- argumentos_extras
    if (is_internal_call) {
      argumentos_internos$internal_call <- NULL
    }

    funcao_objetivo <- function(parametros, ...) {
      res <- tryCatch(-do.call(log_posterior, c(list(parametros), argumentos_internos)), error = function(e) Inf)
      if (length(res) != 1) return(Inf)
      return(res)
    }

    funcao_suporte <- function(parametros, ...) {
      res <- tryCatch(do.call(log_posterior, c(list(parametros), argumentos_internos)), error = function(e) -Inf)
      return(all(is.finite(res)))
    }

    if (!funcao_suporte(x0) || !funcao_suporte(xp0)) {
      stop("Pontos iniciais fora do suporte (log-posterior e -Inf ou gera erro).")
    }

    valor_U <- funcao_objetivo(x0); valor_Up <- funcao_objetivo(xp0)
    x_atual <- x0; xp_atual <- xp0

    matriz_amostras_x <- matrix(NA, nrow = n_iter, ncol = dimensao)
    matriz_amostras_xp <- matrix(NA, nrow = n_iter, ncol = dimensao)
    total_aceitos <- 0

    use_progress_bar <- !is_internal_call
    if (use_progress_bar) {
      barra_progresso <- utils::txtProgressBar(min = 0, max = n_iter, style = 3, width = 50, char = "=")
    }

    for (iteracao in 1:n_iter) {
      movimento <- do.call(executar_movimento, c(
        list(dimensao = dimensao, funcao_objetivo = funcao_objetivo, funcao_suporte = funcao_suporte,
             x = x_atual, valor_U = valor_U, x_prime = xp_atual, valor_Up = valor_Up),
        argumentos_internos
      ))

      if (stats::runif(1) < movimento$alpha) {
        x_atual <- movimento$y; valor_U <- movimento$proposta_U
        xp_atual <- movimento$y_prime; valor_Up <- movimento$proposta_Up
        total_aceitos <- total_aceitos + 1
      }

      matriz_amostras_x[iteracao,] <- x_atual
      matriz_amostras_xp[iteracao,] <- xp_atual
      if (use_progress_bar) {
        utils::setTxtProgressBar(barra_progresso, iteracao)
      }
    }
    if (use_progress_bar) {
      close(barra_progresso)
    }

    taxa_aceitacao <- total_aceitos / n_iter
    if (use_progress_bar) {
      cat(sprintf("\nTaxa de aceitacao: %.2f%%\n", taxa_aceitacao * 100))
    }

    return(list(
      todas_amostras = rbind(matriz_amostras_x, matriz_amostras_xp),
      taxa_aceitacao = taxa_aceitacao,
      num_iteracoes = n_iter,
      dimensao = dimensao
    ))
  }

  # --- BLOCO PARALELO ---
  else {

    if (is.null(n_cores)) {
      n_cores <- max(1, parallel::detectCores() - 1)
    }
    n_cores_usados <- min(n_cadeias, n_cores)

    cat(sprintf("--- Executando T-Walk em modo PARALELO (%d cadeias em %d nucleos) ---\n", n_cadeias, n_cores_usados))

    cl <- parallel::makeCluster(n_cores_usados)
    on.exit(parallel::stopCluster(cl))

    all_objects <- ls(.GlobalEnv)
    parallel::clusterExport(cl, varlist = all_objects, envir = .GlobalEnv)

    parallel::clusterEvalQ(cl, {
      library(mvtnorm)
    })

    executar_uma_cadeia_paralela <- function(indice_cadeia) {
      set.seed(Sys.time() + indice_cadeia)

      dimensao <- length(x0)
      x0_i <- stats::rnorm(dimensao, mean = x0, sd = 0.1)
      xp0_i <- stats::rnorm(dimensao, mean = xp0, sd = 0.1)

      # Usa 'do.call' para construir a chamada de forma segura,
      # passando os argumentos extras (...) corretamente.
      resultado_cadeia <- do.call(twalk, c(
        list(log_posterior = log_posterior, n_iter = n_iter, x0 = x0_i, xp0 = xp0_i,
             n_cadeias = 1, internal_call = TRUE),
        argumentos_extras
      ))
      return(resultado_cadeia)
    }

    cat("Distribuindo trabalho entre os nucleos...\n")
    lista_resultados <- parallel::parLapply(cl, 1:n_cadeias, executar_uma_cadeia_paralela)

    cat("Cadeias concluidas. Combinando resultados...\n")

    amostras_combinadas <- do.call(rbind, lapply(lista_resultados, function(res) res$todas_amostras))
    taxa_aceitacao_media <- mean(sapply(lista_resultados, function(res) res$taxa_aceitacao))
    cat(sprintf("\nTaxa de aceitacao media entre as cadeias: %.2f%%\n", taxa_aceitacao_media * 100))

    return(list(
      todas_amostras = amostras_combinadas,
      taxa_aceitacao = taxa_aceitacao_media,
      num_iteracoes_total = n_iter * n_cadeias,
      dimensao = length(x0),
      cadeias_individuais = lista_resultados
    ))
  }
}
