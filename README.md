# twalkR: An MCMC Sampler Using the t-walk Algorithm

`twalkR` é uma implementação em R do **t-walk**, um amostrador MCMC (Markov Chain Monte Carlo) de propósito geral para distribuições contínuas, ideal para problemas de inferência Bayesiana.

## Visão Geral

Este pacote fornece uma implementação do algoritmo t-walk, conforme proposto originalmente por Christen & Fox (2010). O t-walk é um amostrador MCMC robusto e autoajustável, o que significa que ele não requer um ajuste manual demorado de parâmetros de proposta. Ele é projetado para explorar eficientemente uma ampla gama de distribuições-alvo, mantendo um bom desempenho mesmo em problemas de alta dimensionalidade ou com múltiplas modas.

## Instalação

Tu podes instalar `twalkR` diretamente do GitHub com o pacote `devtools`:

```r
# Se tu não tiveres o pacote 'devtools', instale-o primeiro
# install.packages("devtools")

devtools::install_github("rodrigosqrt3/twalkR", build_vignettes = TRUE)
```

## Exemplo de Uso

Aqui está um exemplo simples de como usar o `twalkR` para amostrar de uma distribuição bimodal.

```r
library(twalkR)

# Definir a log-densidade posterior da distribuição-alvo, neste caso, 
# uma mistura de duas normais
log_posterior_bimodal <- function(x) {
  log(0.5 * dnorm(x, mean = -3, sd = 0.5) + 0.5 * dnorm(x, mean = 3, sd = 0.5))
}

ponto_inicial_1 <- -3
ponto_inicial_2 <- 3

resultado <- twalk(
  log_posterior = log_posterior_bimodal,
  n_iter = 50000,
  x0 = ponto_inicial_1,
  xp0 = ponto_inicial_2
)

burnin <- nrow(resultado$todas_amostras) * 0.2
amostras <- resultado$todas_amostras[-(1:burnin), ]

par(mfrow = c(1, 2))
hist(amostras, breaks = 50, freq = FALSE, 
     main = "Distribuição Posterior", xlab = "Valor do Parâmetro")
lines(density(amostras), col = "blue", lwd = 2)

plot(amostras, type = 'l', col = "grey30", 
     main = "Trace Plot", xlab = "Iteração")
```

## Citação

Este pacote é uma implementação do algoritmo descrito no seguinte artigo. Se você usar `twalkR` em sua pesquisa, por favor, cite o trabalho original:

> Christen, J. A., & Fox, C. (2010). A general purpose sampling algorithm for continuous distributions (the t-walk). *Bayesian Analysis*, 5(2), 263-282. [doi:10.1214/10-BA603](https://doi.org/10.1214/10-BA603)

## Licença

Este pacote está licenciado sob a GPL-3. Veja o arquivo `LICENSE` para mais detalhes.