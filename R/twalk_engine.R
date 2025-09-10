# Funções auxiliares básicas
produto_interno <- function(x) { sum(x*x) }
produto_escalar <- function(x, y) { sum(x*y) }

log_2pi <- log(2*pi)
log_3 <- log(3)

#' Simula o parâmetro beta para o kernel Traverse
#'
#' Amostra beta da distribuição proposta na Seção 2.2 de Christen & Fox (2010).
#' @param at Parâmetro de forma `a_t` (padrão 6).
#' @return Um valor numérico para beta.
#' @keywords internal
simular_beta <- function(at) {
  if (stats::runif(1) < (at - 1) / (2 * at))
    exp(1 / (at + 1) * log(stats::runif(1)))
  else
    exp(1 / (1 - at) * log(stats::runif(1)))
}

#' Gera proposta para o kernel Traverse
#' @param dimensao Dimensão do espaço de parâmetros.
#' @param prob_phi Probabilidade de atualizar cada coordenada.
#' @param x Ponto atual `x`.
#' @param x_prime Ponto atual `x'`.
#' @param beta Parâmetro de passo beta.
#' @return Uma lista com `proposta` e `n_phi` (número de coordenadas movidas).
#' @keywords internal
kernel_traverse <- function(dimensao, prob_phi, x, x_prime, beta) {
  phi <- (stats::runif(dimensao) < prob_phi)
  proposta <- x_prime + beta * (x_prime - x)
  proposta[!phi] <- x[!phi]
  list(proposta = proposta, n_phi = sum(phi))
}

#' Gera proposta para o kernel Walk
#' @param dimensao Dimensão do espaço de parâmetros.
#' @param prob_phi Probabilidade de atualizar cada coordenada.
#' @param parametro_aw Parâmetro de escala `a_w` (padrão 1.5).
#' @param x Ponto atual `x`.
#' @param x_prime Ponto atual `x'`.
#' @return Uma lista com `proposta` e `n_phi`.
#' @keywords internal
kernel_walk <- function(dimensao, prob_phi, parametro_aw, x, x_prime) {
  u <- stats::runif(dimensao)
  phi <- (stats::runif(dimensao) < prob_phi)
  z <- (parametro_aw / (1 + parametro_aw)) * (parametro_aw * u^2 + 2 * u - 1)
  z <- z * phi
  list(proposta = x + (x - x_prime) * z, n_phi = sum(phi))
}

#' Gera proposta para o kernel Blow
#' @param dimensao Dimensão do espaço de parâmetros.
#' @param prob_phi Probabilidade de atualizar cada coordenada.
#' @param x Ponto atual `x`.
#' @param x_prime Ponto atual `x'`.
#' @return Uma lista com `proposta`, `n_phi` e `phi`.
#' @keywords internal
kernel_blow <- function(dimensao, prob_phi, x, x_prime) {
  phi <- (stats::runif(dimensao) < prob_phi)
  sigma <- max(phi * abs(x_prime - x))
  if (sigma < .Machine$double.eps) sigma <- .Machine$double.eps # Evitar sigma = 0

  ruido <- stats::rnorm(dimensao) * phi
  proposta <- x_prime * phi + sigma * ruido + x * (1 - phi)

  list(proposta = proposta, n_phi = sum(phi), phi = phi)
}

#' Calcula o log da densidade da proposta para o kernel Blow
#' @param n_phi Número de coordenadas movidas.
#' @param phi Vetor booleano indicando as coordenadas movidas.
#' @param h Ponto proposto.
#' @param x Ponto de referência.
#' @param x_prime Ponto de referência.
#' @return O valor de -log(g(h)).
#' @keywords internal
log_densidade_blow <- function(n_phi, phi, h, x, x_prime) {
  sigma <- max(phi * abs(x_prime - x))
  if (sigma < .Machine$double.eps) sigma <- .Machine$double.eps

  if (n_phi > 0) {
    diferenca <- h - x_prime
    (n_phi / 2) * log(2 * pi) + n_phi * log(sigma) + 0.5 * sum(diferenca^2) / (sigma^2)
  } else {
    0
  }
}

#' Gera proposta para o kernel Hop
#' @param dimensao Dimensão do espaço de parâmetros.
#' @param prob_phi Probabilidade de atualizar cada coordenada.
#' @param x Ponto atual `x`.
#' @param x_prime Ponto atual `x'`.
#' @return Uma lista com `proposta`, `n_phi` e `phi`.
#' @keywords internal
kernel_hop <- function(dimensao, prob_phi, x, x_prime) {
  phi <- (stats::runif(dimensao) < prob_phi)
  sigma <- max(phi * abs(x_prime - x)) / 3
  if (sigma < .Machine$double.eps) sigma <- .Machine$double.eps

  ruido <- stats::rnorm(dimensao) * phi
  proposta <- x + sigma * ruido
  proposta[!phi] <- x[!phi]

  list(proposta = proposta, n_phi = sum(phi), phi = phi)
}

#' Calcula o log da densidade da proposta para o kernel Hop
#' @param n_phi Número de coordenadas movidas.
#' @param phi Vetor booleano indicando as coordenadas movidas.
#' @param h Ponto proposto.
#' @param x Ponto de referência.
#' @param x_prime Ponto de referência.
#' @return O valor de -log(g(h)).
#' @keywords internal
log_densidade_hop <- function(n_phi, phi, h, x, x_prime) {
  sigma <- max(phi * abs(x_prime - x)) / 3
  if (sigma < .Machine$double.eps) sigma <- .Machine$double.eps

  if (n_phi > 0) {
    diferenca <- h - x
    (n_phi / 2) * log_2pi - n_phi * log_3 + n_phi * log(sigma) + 0.5 * sum(diferenca^2) / (sigma^2)
  } else {
    0
  }
}

#' Executa um único movimento (passo) do algoritmo t-walk
#'
#' Seleciona um dos quatro kernels, gera uma proposta e calcula
#' a probabilidade de aceitação de Metropolis-Hastings.
#' @param ... Argumentos passados para `funcao_objetivo` e `funcao_suporte`.
#' @return Uma lista contendo a proposta e a probabilidade de aceitação.
#' @keywords internal
executar_movimento <- function(dimensao, funcao_objetivo, funcao_suporte, x, valor_U, x_prime, valor_Up,
                               parametro_at = 6, parametro_aw = 1.5,
                               prob_phi = min(dimensao, 4) / dimensao,
                               freq_1 = 0.4918, freq_2 = 0.9836, freq_3 = 0.9918, ...) {

  # Seleção do kernel
  selecao_kernel <- stats::runif(1)

  if (selecao_kernel < freq_1) {	# Kernel Traverse
    direcao <- stats::runif(1)
    if (direcao < 0.5) {
      beta <- simular_beta(parametro_at)
      resultado <- kernel_traverse(dimensao, prob_phi, x_prime, x, beta)
      y_prime <- resultado$proposta; n_phi <- resultado$n_phi; y <- x; proposta_U <- valor_U
      if (funcao_suporte(y_prime, ...)) {
        proposta_Up <- funcao_objetivo(y_prime, ...); if (n_phi == 0) alpha <- 1 else alpha <- exp((valor_U - proposta_U) + (valor_Up - proposta_Up) + (n_phi-2)*log(beta))
      } else { proposta_Up <- NULL; alpha <- 0 }
    } else {
      beta <- simular_beta(parametro_at)
      resultado <- kernel_traverse(dimensao, prob_phi, x, x_prime, beta)
      y <- resultado$proposta; n_phi <- resultado$n_phi; y_prime <- x_prime; proposta_Up <- valor_Up
      if (funcao_suporte(y, ...)) {
        proposta_U <- funcao_objetivo(y, ...); if (n_phi == 0) alpha <- 1 else alpha <- exp((valor_U - proposta_U) + (valor_Up - proposta_Up) + (n_phi-2)*log(beta))
      } else { proposta_U <- NULL; alpha <- 0 }
    }
  } else if (selecao_kernel < freq_2) { # Kernel Walk
    direcao <- stats::runif(1)
    if (direcao < 0.5) {
      resultado <- kernel_walk(dimensao, prob_phi, parametro_aw, x_prime, x)
      y_prime <- resultado$proposta; n_phi <- resultado$n_phi; y <- x; proposta_U <- valor_U
      if (funcao_suporte(y_prime, ...) && (all(abs(y_prime - y) > 0))) {
        proposta_Up <- funcao_objetivo(y_prime, ...); alpha <- exp((valor_U - proposta_U) + (valor_Up - proposta_Up))
      } else { proposta_Up <- NULL; alpha <- 0 }
    } else {
      resultado <- kernel_walk(dimensao, prob_phi, parametro_aw, x, x_prime)
      y <- resultado$proposta; n_phi <- resultado$n_phi; y_prime <- x_prime; proposta_Up <- valor_Up
      if (funcao_suporte(y, ...) && (all(abs(y_prime - y) > 0))) {
        proposta_U <- funcao_objetivo(y, ...); alpha <- exp((valor_U - proposta_U) + (valor_Up - proposta_Up))
      } else { proposta_U <- NULL; alpha <- 0 }
    }
  } else if (selecao_kernel < freq_3) { # Kernel Blow
    direcao <- stats::runif(1)
    if (direcao < 0.5) {
      resultado <- kernel_blow(dimensao, prob_phi, x_prime, x)
      y_prime <- resultado$proposta; n_phi <- resultado$n_phi; phi <- resultado$phi; y <- x; proposta_U <- valor_U
      if (funcao_suporte(y_prime, ...) && all(y_prime != x)) {
        proposta_Up <- funcao_objetivo(y_prime, ...); W1 <- log_densidade_blow(n_phi, phi, y_prime, x_prime, x); W2 <- log_densidade_blow(n_phi, phi, x_prime, y_prime, x); alpha <- exp((valor_U - proposta_U) + (valor_Up - proposta_Up) + (W1 - W2))
      } else { proposta_Up <- NULL; alpha <- 0 }
    } else {
      resultado <- kernel_blow(dimensao, prob_phi, x, x_prime)
      y <- resultado$proposta; n_phi <- resultado$n_phi; phi <- resultado$phi; y_prime <- x_prime; proposta_Up <- valor_Up
      if (funcao_suporte(y, ...) && all(y != x_prime)) {
        proposta_U <- funcao_objetivo(y, ...); W1 <- log_densidade_blow(n_phi, phi, y, x, x_prime); W2 <- log_densidade_blow(n_phi, phi, x, y, x_prime); alpha <- exp((valor_U - proposta_U) + (valor_Up - proposta_Up) + (W1 - W2))
      } else { proposta_U <- NULL; alpha <- 0 }
    }
  } else { # Kernel Hop
    direcao <- stats::runif(1)
    if (direcao < 0.5) {
      resultado <- kernel_hop(dimensao, prob_phi, x_prime, x)
      y_prime <- resultado$proposta; n_phi <- resultado$n_phi; phi <- resultado$phi; y <- x; proposta_U <- valor_U
      if (funcao_suporte(y_prime, ...) && all(y_prime != x)) {
        proposta_Up <- funcao_objetivo(y_prime, ...); W1 <- log_densidade_hop(n_phi, phi, y_prime, x_prime, x); W2 <- log_densidade_hop(n_phi, phi, x_prime, y_prime, x); alpha <- exp((valor_U - proposta_U) + (valor_Up - proposta_Up) + (W1 - W2))
      } else { proposta_Up <- NULL; alpha <- 0 }
    } else {
      resultado <- kernel_hop(dimensao, prob_phi, x, x_prime)
      y <- resultado$proposta; n_phi <- resultado$n_phi; phi <- resultado$phi; y_prime <- x_prime; proposta_Up <- valor_Up
      if (funcao_suporte(y, ...) && all(y != x_prime)) {
        proposta_U <- funcao_objetivo(y, ...); W1 <- log_densidade_hop(n_phi, phi, y, x, x_prime); W2 <- log_densidade_hop(n_phi, phi, x, y, x_prime); alpha <- exp((valor_U - proposta_U) + (valor_Up - proposta_Up) + (W1 - W2))
      } else { proposta_U <- NULL; alpha <- 0 }
    }
  }

  if (is.nan(alpha)) { alpha <- 0 }

  list(y = y, proposta_U = proposta_U, y_prime = y_prime, proposta_Up = proposta_Up, alpha = alpha)
}
