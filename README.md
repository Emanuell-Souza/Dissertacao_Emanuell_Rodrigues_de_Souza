# Predição Clínica Pan-Câncer e Vulnerabilidades Terapêuticas em Vias de Morte Celular Regulada

<p align="center">
  <img src="https://img.shields.io/badge/Linguagem-R-276DC3?style=for-the-badge&logo=r" alt="R">
  <img src="https://img.shields.io/badge/Domínio-Bioinformática-2E8B57?style=for-the-badge" alt="Bioinformatics">
  <img src="https://img.shields.io/badge/IA-SuperLearner%20%2B%20SHAP-FF6F00?style=for-the-badge" alt="AI">
  <img src="https://img.shields.io/badge/Status-Dissertação%20em%20Andamento-yellow?style=for-the-badge" alt="Status">
</p>

---

## 🎓 Sobre a Dissertação

Este repositório contém todos os scripts e rotinas analíticas desenvolvidos para a dissertação de mestrado:

> **Emanuell Rodrigues de Souza**  
> *Predição Clínica Pan-Câncer e Vulnerabilidades Terapêuticas em Vias de Morte Celular Regulada*  
> Dissertação apresentada ao Centro de Biociências e Biotecnologia da **Universidade Estadual do Norte Fluminense Darcy Ribeiro (UENF)** como parte das exigências para obtenção do título de **Mestre em Biociências e Biotecnologia**.  
> **Orientador:** Prof. Dr. Enrique Medina-Acosta  
> **Programa:** Pós-Graduação em Biociências e Biotecnologia (PGBB) – Laboratório de Biotecnologia (LBT)  
> **Local:** Campos dos Goytacazes – RJ, 2026

---

## 📄 Resumo

O câncer representa uma das principais causas de mortalidade no mundo e caracteriza-se por elevada heterogeneidade molecular, celular e clínica, dificultando a identificação de biomarcadores prognósticos e de vulnerabilidades terapêuticas aplicáveis à medicina de precisão. Nesse contexto, **vias de morte celular regulada (RCD)** desempenham papel central na progressão tumoral, evasão terapêutica e remodelação do microambiente tumoral.

Esta dissertação teve como objetivo desenvolver uma abordagem **Pan-câncer** baseada em **integração multiômica** e **aprendizado de máquina** para identificar assinaturas associadas à predição clínica e a potenciais vulnerabilidades terapêuticas relacionadas à RCD. Foram integrados dados transcriptômicos e genômicos provenientes das bases **TCGA**, **CCLE** e **DepMap**.

Os resultados demonstraram desempenho preditivo consistente, além da identificação de assinaturas associadas a compartimentos celulares específicos, dependências funcionais críticas e perfis distintos de sensibilidade farmacológica. Em conjunto, os achados reforçam o potencial da integração entre multiômica, aprendizado de máquina e inteligência artificial explicável para a descoberta de biomarcadores prognósticos e alvos terapêuticos em oncologia.

**Palavras-chave:** aprendizado de máquina; câncer; células únicas; essencialidade gênica; morte celular regulada; sensibilidade a drogas; transcriptômica espacial.

---

## 📂 Estrutura do Repositório

O repositório está organizado em módulos que refletem as etapas sequenciais do pipeline analítico:

| Pasta | Conteúdo |
|---|---|
| `1-Harmonização_imputação/` | Harmonização de variáveis clínicas e moleculares; estratégias de imputação de dados faltantes |
| `2-RSF/` | Modelos de sobrevivência baseados em **Random Survival Forest (RSF)** |
| `3-XGB/` | Modelos de sobrevivência com **XGBoost** |
| `4-superlearner/` | Ensemble **SuperLearner** combinando RSF, XGBoost e outros estimadores base |
| `5-scRNA/` | Análise de **Transcriptômica de Célula Única (scRNA-seq)** |
| `6-stRNA/` | Análise de **Transcriptômica Espacial (Spatial Transcriptomics)** |
| `7-Essencialidade/` | Avaliação de **essencialidade gênica** via dados CRISPR/RNAi do DepMap |
| `8-Sensibilidade/` | Integração com dados de **sensibilidade a drogas** (PRISM – Broad Institute) |

---

## 🔬 Pipeline Analítico

```
Dados TCGA / CCLE / DepMap
         │
         ▼
 1. Harmonização & Imputação
         │
         ▼
 2. Modelos de Sobrevivência
    ├── Random Survival Forest (RSF)
    ├── XGBoost
    └── SuperLearner (Ensemble)
         │
         ▼
 3. IA Explicável (SHAP / LIME)
    └── Identificação de Assinaturas Moleculares
         │
         ├──► 4. scRNA-seq (distribuição celular)
         ├──► 5. Transcriptômica Espacial (distribuição espacial)
         ├──► 6. Essencialidade Gênica (DepMap/CRISPR)
         └──► 7. Sensibilidade a Drogas (PRISM)
```

---

## 🎯 Desfechos Clínicos Avaliados

- **OS** – Overall Survival (Sobrevida Global)
- **DSS** – Disease-Specific Survival (Sobrevida Doença-Específica)
- **DFI** – Disease-Free Interval (Intervalo Livre de Doença)
- **PFI** – Progression-Free Interval (Intervalo Livre de Progressão)

---

## 🛠️ Tecnologias Utilizadas

- **Linguagem:** R
- **Modelos:** `randomForestSRC`, `xgboost`, `SuperLearner`
- **XAI:** `SHAPforxgboost`, `lime`
- **Visualização:** `ggplot2`, `Seurat`, `ggpubr`
- **Bases de Dados:** TCGA, CCLE, DepMap (Broad Institute)

---

## 📜 Como Citar

Se você utilizar este código ou os resultados em sua pesquisa, por favor cite:

> Rodrigues de Souza, E. **Predição Clínica Pan-Câncer e Vulnerabilidades Terapêuticas em Vias de Morte Celular Regulada.** Dissertação (Mestrado em Biociências e Biotecnologia) – Universidade Estadual do Norte Fluminense Darcy Ribeiro, Campos dos Goytacazes, 2026.

---

## 📧 Contato

- **Autor:** Emanuell Rodrigues de Souza
- **E-mail:** souza.ers00@gmail.com
- **Orientador:** Prof. Dr. Enrique Medina-Acosta – UENF/LBT

---


## 💰 Financiamento

Este trabalho recebeu apoio institucional do **Programa de Apoio à Pesquisa Institucional (PAPIC; Processo UENF001/2024)** e do **PAPIC PLUS (Processo UENF001/2025)**. O financiamento adicional de infraestrutura foi fornecido pela **Financiadora de Estudos e Projetos (Finep)** e pelo **Fundo Nacional de Desenvolvimento Científico e Tecnológico (FNDCT)** no âmbito do programa **PROINFRA 2021** (Convênio nº 0.1.22.0442.00).

Emanuell Rodrigues de Souza foi contemplado com **bolsa de mestrado da Coordenação de Aperfeiçoamento de Pessoal de Nível Superior (CAPES), Brasil**.
