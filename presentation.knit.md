---
title: "Simulation Based Methods for Network Inference"
author: "Marthyna Luiza WEBER"
institute: "Grenoble INP - Ensimag"
date: "2023-05-15"
output:
  beamer_presentation
---
\tableofcontents

# Initial definitions

## 1. Partial Correlation

Correlation coefficient between $X_1$ and $X_2$ after removing the influence of $Y$, accounting for the scaling of the variables.

-   Partial covariance: The partial covariance between $X_1$ and $X_2$ with reference to $Y$ is calculated with:

$$
Cov(X_1, X_2 \cdot Y) = \mathbb{E}[(X_1 - \hat{X}1(Y))(X_2 - \hat{X}2(Y))]
$$

-   $\hat{X}_l(Y)$ is the projection of $X_l$ on $Y$, which is the expected value of $X_l$ given $Y$:

$$
\hat{X}_l(Y) = \mathbb{E}(X_l)+\frac{Cov(X_l,Y)}{Var(Y)}(Y-\mathbb{E}(Y))
$$

## 1. Partial correlation

-   Lemma: $$ Cov(X_\mathit{l},X_\mathit{m}\cdot Y) = Cov(X_\mathit{l},X_\mathit{m}) - Cov(X_\mathit{l},Y)Var(Y)^{-1}Cov(X_\mathit{m}Y)^T $$ \small Given a vector $\mathbf{X} = (X_1, \ldots, X_p)$ of $p$ random variables. We have:

$$ Cov^{\text{partial}}(\mathbf{X}_{\mathit{l,m}}) = Cov(X_\mathit{l},X_\mathit{m}\cdot X_{\mathit{V\backslash\{l,m\}}}) $$

-   Partial correlation: $$ \rho^{\text{partial}}(\mathbf{X}_{\mathit{l,m}})=\frac{Cov(X_\mathit{l},X_\mathit{m}\cdot X_{\mathit{V\backslash\{l,m\}}})}{\sqrt{Var(X_\mathit{l} \cdot X_{\mathit{V\backslash\{l,m\}}}) Var(X_\mathit{m} \cdot X_{\mathit{V\backslash\{l,m\}}})}} $$

## 2. Regular vines

\small A vine in which two edges in tree $T_i$ are joined by an edge in tree $T_{i+1}$ only if these edges share a common node.

\begin{columns}[T]
  \begin{column}{.5\textwidth}
    \centering
    \includegraphics[width=1.1\textwidth]{regular_vine.png}
  \end{column}
  \begin{column}{.5\textwidth}
    \begin{itemize}
      \item Constraint set $U_e$: variables reachable from edge $e$
      \item Conditioning set $D_e$: variables shared between $U_e$ and edge adjacent to $e$ in the next tree
      \item Conditioned set $\{C_{1e}, C_{2e}\}$: symmetric difference of $U_e$ and $D_e$
    \end{itemize}
  \end{column}
\end{columns}

- $\{L|K\}$: the constraint set, with conditioned set $L$ and conditioning set $K$.

## 4. C-vines

- \small C-vines: vines where each tree $T_i$ has a unique node of degree $d-i$
- \small D-vines: vines where each node in $T_i$ has a degree at most 2.

\begin{columns}[T]
  \begin{column}{0.5\textwidth}
    \centering
    \includegraphics[width=\textwidth]{c-vines.png}
    \tiny A C-vine
  \end{column}
  \begin{column}{0.5\textwidth}
    \centering
    \includegraphics[width=\textwidth]{d-vines.png}
    \tiny A D-vine
  \end{column}
\end{columns}

## 5. Generatating random correlation matrices with C-vines
Algorithm to generate a random correlation matrix $\boldsymbol{R}$ with density proportional to $det(\boldsymbol{R})^{\eta-1}$, with $\eta > 1$:
\begin{enumerate}
    \item Initialize $\beta = \eta + \frac{d-1}{2}$
    \item Loop for $k = 1, \ldots, d-1$:
    \begin{enumerate}
        \item $\beta = \beta - \frac{1}{2}$
        \item Loop for $i = k+1, \ldots, d$:
        \begin{enumerate}
            \item generate $p_{k,i;1,\ldots,k-1}$ $\sim$ Beta$(\beta, \beta)$  on $(-1,1)$
            \item use the recursive formula for partial correlations calculation
        \end{enumerate}
    \end{enumerate}
    \item Return $\boldsymbol{R}$, a $d \times d$ correlation matrix
\end{enumerate}

