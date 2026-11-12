
# ============================================================
# MATA DO MAMÃO — REANÁLISE ESTATÍSTICA
# ============================================================
# Entrada: Todos_Com_Indices_Spectrais (1).csv
# Saída: pasta resultados_mata_mamao_R
# ============================================================

# 1. PACOTES
pacotes <- c("ggplot2", "readr", "dplyr", "tidyr", "tibble")

faltam <- pacotes[!sapply(pacotes, requireNamespace, quietly = TRUE)]
if (length(faltam) > 0) install.packages(faltam)

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(tibble)

# 2. CONFIGURAÇÕES
# Ao executar, selecione o CSV na janela de arquivos.
arquivo <- file.choose()

B <- 5000L
set.seed(20261112)

saida <- "resultados_mata_mamao_R"
dir.create(saida, showWarnings = FALSE, recursive = TRUE)

brutos <- read_csv(arquivo, show_col_types = FALSE)

colunas <- c(
  "geom_tipo", "id_transecto", "total_pnt",
  "p_lenhoso", "p_nao_len", "p_solo_nu",
  "p_len_nlen", "p_nlen_len",
  paste0(c("NDVI", "GNDVI", "SAVI",
           "BSI", "NDMI", "NBR"), "_mean")
)

faltantes <- setdiff(colunas, names(brutos))
if (length(faltantes)) {
  stop("Colunas ausentes: ", paste(faltantes, collapse = ", "))
}

geoms <- c(
  "Pontos", "Linhas_50m", "Buffers_25m",
  "Pixels_S2_10m", "Pixels_S2_20m", "Pixels_L8_30m"
)

rot_geoms <- c(
  "Pontos agregados", "Linhas 50 m", "Buffers 25 m",
  "Sentinel-2 10 m", "Sentinel-2 20 m", "Landsat 8 30 m"
)

classes <- c("p_lenhoso", "p_nao_len", "p_solo_nu", "p_dupla")

rot_classes <- c(
  "Lenhoso isolado", "Não lenhoso isolado",
  "Solo nu", "Classe dupla"
)

indices <- paste0(
  c("NDVI", "GNDVI", "SAVI", "BSI", "NDMI", "NBR"),
  "_mean"
)

rot_indices <- sub("_mean$", "", indices)

# ============================================================
# 3. AUDITORIA DOS DADOS
# ============================================================

auditoria <- brutos %>%
  group_by(geom_tipo) %>%
  summarise(
    unidades = n(),
    transectos = n_distinct(id_transecto),
    sem_pontos = sum(is.na(total_pnt) | total_pnt == 0),
    mediana_pontos = median(total_pnt, na.rm = TRUE),
    soma_pontos = sum(total_pnt, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(auditoria, file.path(saida, "01_auditoria.csv"))

print(auditoria)

# Excluir unidades sem registros de campo
dados <- brutos %>%
  mutate(p_dupla = p_len_nlen + p_nlen_len) %>%
  filter(
    !is.na(total_pnt),
    total_pnt > 0,
    !is.na(id_transecto),
    geom_tipo %in% geoms
  )

# ============================================================
# 4. CORRELAÇÕES DE PEARSON
# ============================================================

cor_segura <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  
  if (sum(ok) < 3L ||
      sd(x[ok]) == 0 ||
      sd(y[ok]) == 0) {
    return(NA_real_)
  }
  
  unname(cor(x[ok], y[ok], method = "pearson"))
}

cor_p <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  
  if (sum(ok) < 3L ||
      sd(x[ok]) == 0 ||
      sd(y[ok]) == 0) {
    return(NA_real_)
  }
  
  cor.test(x[ok], y[ok], method = "pearson")$p.value
}

comb <- expand.grid(
  geom_tipo = geoms,
  classe = classes,
  indice = indices,
  stringsAsFactors = FALSE
)

correlacoes <- bind_rows(
  lapply(seq_len(nrow(comb)), function(i) {
    z <- comb[i, ]
    d <- dados[dados$geom_tipo == z$geom_tipo, , drop = FALSE]
    
    x <- d[[z$classe]]
    y <- d[[z$indice]]
    
    tibble(
      geom_tipo = z$geom_tipo,
      classe = z$classe,
      indice = z$indice,
      n = sum(is.finite(x) & is.finite(y)),
      r = cor_segura(x, y),
      p_pearson_descritivo = cor_p(x, y)
    )
  })
)

write_csv(
  correlacoes,
  file.path(saida, "02_todas_correlacoes.csv")
)

# ============================================================
# 5. MATRIZ DE CORRELAÇÕES — FIGURA DA NOTA
# ============================================================

matriz <- correlacoes %>%
  mutate(
    geometria = factor(
      geom_tipo, levels = geoms, labels = rot_geoms
    ),
    classe_nome = factor(
      classe,
      levels = rev(classes),
      labels = rev(rot_classes)
    ),
    indice_nome = factor(
      sub("_mean$", "", indice),
      levels = rot_indices
    )
  )

grafico <- ggplot(
  matriz,
  aes(x = indice_nome, y = classe_nome, fill = r)
) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(
    aes(label = ifelse(is.na(r), "NA", sprintf("%.2f", r))),
    size = 3.2
  ) +
  facet_wrap(~geometria, ncol = 2) +
  scale_fill_gradient2(
    low = "#B44B45",
    mid = "white",
    high = "#28734B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson (r)"
  ) +
  labs(
    title = "Matriz de correlações entre cobertura vegetal e índices espectrais",
    subtitle = "Mata do Mamão | unidades sem registros de campo excluídas",
    x = "Índice espectral",
    y = "Classe de cobertura"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

print(grafico)

ggsave(
  file.path(saida, "03_matriz_correlacao.png"),
  grafico, width = 12, height = 11, dpi = 300, bg = "white"
)

ggsave(
  file.path(saida, "03_matriz_correlacao.pdf"),
  grafico, width = 12, height = 11, bg = "white"
)

# ============================================================
# 6. TABELA 1 DA NOTA CIENTÍFICA
# ============================================================

classes_tabela <- c("p_lenhoso", "p_dupla", "p_solo_nu")

tabela_longa <- correlacoes %>%
  filter(classe %in% classes_tabela, !is.na(r)) %>%
  group_by(geom_tipo, classe) %>%
  arrange(desc(abs(r)), indice, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    geometria = factor(
      geom_tipo, levels = geoms, labels = rot_geoms
    ),
    classe_nome = factor(
      classe,
      levels = classes_tabela,
      labels = c("Lenhoso", "Classe dupla", "Solo nu")
    ),
    valor = sprintf(
      "%s (%.3f)", sub("_mean$", "", indice), r
    )
  ) %>%
  select(geometria, classe_nome, valor, r, indice, n)

write_csv(
  tabela_longa,
  file.path(saida, "04_tabela_nota_detalhada.csv")
)

tabela_nota <- tabela_longa %>%
  select(geometria, classe_nome, valor) %>%
  pivot_wider(
    names_from = classe_nome,
    values_from = valor
  ) %>%
  arrange(geometria)

write_csv(
  tabela_nota,
  file.path(saida, "04_tabela_nota.csv")
)

print(tabela_nota)

# ============================================================
# 7. BOOTSTRAP POR TRANSECTO
# ============================================================

# Três relações centrais da nota científica
pares <- tibble(
  classe = c("p_lenhoso", "p_dupla", "p_solo_nu"),
  indice = c("GNDVI_mean", "GNDVI_mean", "NBR_mean")
)

ids <- unique(dados$id_transecto)
N <- length(ids)

# Organizar os registros por transecto
base_pares <- vector("list", nrow(pares))

for (k in seq_len(nrow(pares))) {
  z <- pares[k, ]
  
  base_pares[[k]] <- lapply(geoms, function(g) {
    d <- dados[
      dados$geom_tipo == g,
      c("id_transecto", z$classe, z$indice)
    ]
    
    list(
      x = d[[z$classe]],
      y = d[[z$indice]],
      por_id = split(
        seq_len(nrow(d)),
        as.character(d$id_transecto)
      )
    )
  })
}

boot <- array(
  NA_real_,
  dim = c(B, length(geoms), nrow(pares)),
  dimnames = list(
    NULL,
    geoms,
    paste(pares$classe, pares$indice, sep = "__")
  )
)

# Reamostrar transectos completos
for (b in seq_len(B)) {
  amostra <- sample(ids, N, replace = TRUE)
  
  for (k in seq_len(nrow(pares))) {
    for (g in seq_along(geoms)) {
      d <- base_pares[[k]][[g]]
      
      pos <- unlist(
        d$por_id[as.character(amostra)],
        use.names = FALSE
      )
      
      boot[b, g, k] <- cor_segura(d$x[pos], d$y[pos])
    }
  }
  
  if (b %% 500 == 0) {
    message("Bootstrap: ", b, " de ", B)
  }
}

ci <- function(v) {
  v <- v[is.finite(v)]
  
  if (length(v) < 100L) {
    return(c(NA_real_, NA_real_))
  }
  
  as.numeric(
    quantile(v, c(0.025, 0.975), names = FALSE)
  )
}

res_boot <- list()
comparacoes <- list()

for (k in seq_len(nrow(pares))) {
  z <- pares[k, ]
  
  # Intervalos de confiança por geometria
  for (g in seq_along(geoms)) {
    d <- dados[dados$geom_tipo == geoms[g], ]
    
    r <- cor_segura(d[[z$classe]], d[[z$indice]])
    intervalo <- ci(boot[, g, k])
    
    res_boot[[length(res_boot) + 1L]] <- tibble(
      classe = z$classe,
      indice = z$indice,
      geometria = geoms[g],
      n = nrow(d),
      r = r,
      ic95_inferior = intervalo[1],
      ic95_superior = intervalo[2]
    )
  }
  
  # Comparar buffer com as outras geometrias
  for (g in which(geoms != "Buffers_25m")) {
    dif <- boot[, "Buffers_25m", k] - boot[, g, k]
    
    intervalo <- ci(dif)
    dif <- dif[is.finite(dif)]
    
    d1 <- dados[dados$geom_tipo == "Buffers_25m", ]
    d2 <- dados[dados$geom_tipo == geoms[g], ]
    
    observado <-
      cor_segura(d1[[z$classe]], d1[[z$indice]]) -
      cor_segura(d2[[z$classe]], d2[[z$indice]])
    
    comparacoes[[length(comparacoes) + 1L]] <- tibble(
      classe = z$classe,
      indice = z$indice,
      comparacao = paste("Buffers -", geoms[g]),
      diferenca_r = observado,
      ic95_inferior = intervalo[1],
      ic95_superior = intervalo[2]
    )
  }
}

res_boot <- bind_rows(res_boot)
comparacoes <- bind_rows(comparacoes)

write_csv(
  res_boot,
  file.path(saida, "05_correlacoes_bootstrap.csv")
)

write_csv(
  comparacoes,
  file.path(saida, "06_contrastes_bootstrap.csv")
)

# ============================================================
# 8. ANÁLISE DE SENSIBILIDADE
# ============================================================

# Número mínimo de pontos de campo por unidade
limiares <- c(1L, 5L, 10L, 15L, 20L, 30L)

sensibilidade <- bind_rows(
  lapply(limiares, function(lim) {
    bind_rows(
      lapply(seq_len(nrow(pares)), function(k) {
        z <- pares[k, ]
        
        bind_rows(
          lapply(geoms, function(g) {
            d <- dados %>%
              filter(
                geom_tipo == g,
                total_pnt >= lim
              )
            
            tibble(
              min_pontos = lim,
              geometria = g,
              classe = z$classe,
              indice = z$indice,
              n = nrow(d),
              transectos = n_distinct(d$id_transecto),
              r = cor_segura(
                d[[z$classe]],
                d[[z$indice]]
              )
            )
          })
        )
      })
    )
  })
)

write_csv(
  sensibilidade,
  file.path(saida, "07_sensibilidade_pontos.csv")
)

# ============================================================
# 9. FINALIZAÇÃO
# ============================================================

nota <- c(
  "MATA DO MAMÃO - REANÁLISE",
  paste("Arquivo:", arquivo),
  paste("Transectos:", N),
  paste("Reamostragens bootstrap:", B),
  "Excluídas unidades com total_pnt=0.",
  "Bootstrap realizado por transecto.",
  "As diferenças entre sensores não representam efeito puro de geometria.",
  "A análise temporal exige nova extração espectral.",
  "Verificar cálculo do SAVI e máscaras de nuvens no GEE."
)

writeLines(
  nota,
  file.path(saida, "08_notas_metodologicas.txt")
)

message("ANÁLISE CONCLUÍDA!")
message("Resultados salvos em: ", normalizePath(saida))

