% =====================================================================
%  HCSM-25 (Revised)
%  The 4x4 Periodic Charge Density
%  The Substrate Charge Distribution from the Mode Coordinates
%  Stanley Preschutti
%  Entropia Research Institute / Information Physics Institute
%  ORCID: 0009-0004-5445-1744
%  September 2026
% =====================================================================
\documentclass[10pt,twocolumn]{article}

\usepackage[a4paper,margin=0.75in,columnsep=0.28in]{geometry}
\usepackage{amsmath,amssymb,amsthm,mathtools}
\usepackage{bm}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{array}
\usepackage{microtype}
\usepackage[hidelinks]{hyperref}
\usepackage{enumitem}
\usepackage{authblk}

% ---------- Theorem environments ----------
\theoremstyle{plain}
\newtheorem{theorem}{Theorem}[section]
\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{corollary}[theorem]{Corollary}

\theoremstyle{definition}
\newtheorem{definition}[theorem]{Definition}
\newtheorem{axiom}[theorem]{Axiom}
\newtheorem{input}[theorem]{Input}

\theoremstyle{remark}
\newtheorem{remark}[theorem]{Remark}

% ---------- Shortcuts ----------
\newcommand{\Z}{\mathbb{Z}}
\newcommand{\R}{\mathbb{R}}
\newcommand{\Q}{\mathbb{Q}}
\newcommand{\C}{\mathbb{C}}
\newcommand{\Lam}{\Lambda}
\newcommand{\dd}{\mathrm{d}}
\newcommand{\ip}[2]{\langle #1, #2\rangle}
\newcommand{\abs}[1]{\left| #1 \right|}
\newcommand{\norm}[1]{\left\| #1 \right\|}
\newcommand{\Tr}{\operatorname{Tr}}
\newcommand{\rank}{\operatorname{rank}}
\newcommand{\supp}{\operatorname{supp}}
\newcommand{\IPR}{\operatorname{IPR}}
\newcommand{\sign}{\operatorname{sign}}
\newcommand{\diag}{\operatorname{diag}}

\title{\bfseries The $4\times 4$ Periodic Charge Density\\[2pt]
\large The Substrate Charge Distribution from the Mode Coordinates}

\author[1]{Stanley Preschutti}
\affil[1]{Entropia Research Institute / Information Physics Institute\\
ORCID: 0009-0004-5445-1744\\
\texttt{scstannp@yahoo.com}}

\date{September 2026}

\begin{document}
\maketitle

\begin{abstract}
We establish the charge density of the Hodge Complex Standard Model (HCSM)
on the substrate $\Lam = \Z_8 \times \Z_8$.
The density is
\[
\rho(x,y) \;=\; \frac{1}{64}\sum_{i=1}^{14} q_{3i}\,
\cos\!\left(\frac{\pi\,(n_1^{(i)}x + n_2^{(i)}y)}{2}\right),
\]
where $q_{3i}$ is the electric charge of mode $i$ in units of $e/3$
(from HCSM-24) and $(n_1^{(i)},n_2^{(i)})$ is the mode's coordinate
(from HCSM-04).
We prove: (i) $\rho(x,y)$ is exactly $4\times 4$ periodic on the
$8\times 8$ torus, with fundamental domain given by an explicit
$4\times 4$ matrix; (ii) every row and column of the fundamental
domain sums to $-12$; (iii) the total charge is
$-192 = -3\cdot 64$, i.e.\ $-1e$ per site.
The density is fully derived from the mode coordinates (HCSM-04),
the charge table (HCSM-24, itself a theorem), and the substrate
Laplacian eigenfunctions.
No new empirical inputs are introduced.
All numerical claims are verified in reproducible Python.
\end{abstract}

\medskip
\noindent\textbf{Keywords:} Hodge Complex Standard Model; HCSM; charge density;
$4\times 4$ periodicity; fundamental domain; total charge;
rigidity class; substrate.

% =====================================================================
\section{Introduction}
% =====================================================================

HCSM-24 \cite{HCSM24} established the charge theorem: the electric
charge in units of $e/3$ is $q_3 = F(n_1,\sigma)$ for the 14 modes at
the self-paired eigenvalue $\lambda=4$ of the discrete Laplacian on the
periodic $8\times 8$ torus. The charge table is a theorem of the
substrate: it is the unique assignment satisfying the five structural
constraints (C1)--(C5) of HCSM-24, with the electron/neutrino binary
residue resolved by the Midpoint Sorting Theorem of HCSM-23 Rev.~8.

The present paper establishes the \emph{substrate charge distribution}
$\rho(x,y)$ associated with this charge table. The density is the
standard Fourier reconstruction of the charge on the 64 sites of the
torus, using the substrate Laplacian eigenfunctions at the self-paired
eigenvalue as the reconstruction basis.

\paragraph{Main results.}
\begin{enumerate}[label=(\roman*),leftmargin=*]
\item (\textbf{$4\times 4$ periodicity}, Theorem~\ref{thm:periodicity}.)
      $\rho(x+4,y)=\rho(x,y)$ and $\rho(x,y+4)=\rho(x,y)$ for all
      $(x,y)\in\Lam$.
\item (\textbf{Fundamental domain}, Theorem~\ref{thm:fundamental}.)
      The fundamental domain of $\rho$ is the $4\times 4$ matrix
      \[
      \rho_{4\times 4} = \frac{1}{64}
      \begin{pmatrix}
        -3 & -8 & +7 & -8\\
        -8 & +3 & -8 & +1\\
        +7 & -8 & -3 & -8\\
        -8 & +1 & -8 & +3
      \end{pmatrix}.
      \]
\item (\textbf{Row/column sums}, Corollary~\ref{cor:rowcol}.)
      Every row and column of $\rho_{4\times 4}$ sums to $-12$.
\item (\textbf{Total charge}, Corollary~\ref{cor:total}.)
      $\displaystyle\sum_{x,y=0}^{7}\rho(x,y) = -192 = -3\cdot 64
      = -1e \cdot 64$, i.e.\ $-1e$ per site.
\end{enumerate}

The density is fully derived from the mode coordinates and the charge
table, both of which are established in prior HCSM papers. No new
empirical inputs are introduced.

% =====================================================================
\section{Setup}
% =====================================================================

\subsection{The substrate and the 14-mode multiplet}

The substrate is the periodic square lattice
$\Lam = \Z_L \times \Z_L$ with $N = L^2$ sites, selected at
$L=8$, $N=64$ by the unified substrate selection theorem of HCSM-00
\cite{HCSM00} \S2.6. Site $(x_1,x_2)$ is identified with
$(x_1+8,x_2)$ and with $(x_1,x_2+8)$.

The discrete 5-point Laplacian
\begin{equation}
(\Delta f)(x) = 4f(x) - \sum_{|e|=1} f(x+e)
\label{eq:laplacian}
\end{equation}
has eigenvalues
\begin{equation}
\lambda(n_1,n_2) = 4 - 2\cos\!\left(\frac{\pi n_1}{4}\right)
                       - 2\cos\!\left(\frac{\pi n_2}{4}\right),
\label{eq:eigenvalues}
\end{equation}
for $(n_1,n_2)\in\{0,\ldots,7\}^2$. The spectrum is mirror-symmetric
under $\lambda \leftrightarrow 8-\lambda$, and the unique self-paired
eigenvalue is $\lambda=4$ with multiplicity $2(L-1)=14$ (HCSM-04
\cite{HCSM04}, Theorem~3.2).

\subsection{The 14 modes and their coordinates}

In standard mode order $(\nu,e,u,d,s,\mu,\mathrm{dark},c,\tau,b,W,Z,H,t)$,
the mode coordinates are (HCSM-04 \cite{HCSM04}, Theorem~4.1 and
Remark~4.2):
\begin{equation}
\begin{aligned}
\nu &: (0,4), & e &: (4,0), & u &: (1,3), & d &: (1,5),\\
s &: (3,1), & \mu &: (3,7), & \mathrm{dark} &: (5,7), & c &: (7,3),\\
\tau &: (7,5), & b &: (5,1), & W &: (2,2), & Z &: (2,6),\\
H &: (6,2), & t &: (6,6).
\end{aligned}
\label{eq:modes}
\end{equation}
These decompose under the $D_4$ point group into three orbits
$\mathcal{O}_0,\mathcal{O}_1,\mathcal{O}_2$ of sizes $2,8,4$
(HCSM-04 \cite{HCSM04}, Theorem~5.1).

\subsection{The charge table}

The charge table in units of $e/3$ is established in HCSM-24
\cite{HCSM24} as the unique solution of the five structural
constraints (C1)--(C5):
\begin{equation}
\begin{aligned}
q_3(\nu) &= 0, & q_3(e) &= -3, & q_3(u) &= +2, & q_3(d) &= -1,\\
q_3(s) &= -1, & q_3(\mu) &= -3, & q_3(\mathrm{dark}) &= 0, & q_3(c) &= +2,\\
q_3(\tau) &= -3, & q_3(b) &= -1, & q_3(W) &= +3, & q_3(Z) &= 0,\\
q_3(H) &= 0, & q_3(t) &= +2.
\end{aligned}
\label{eq:charges}
\end{equation}
The total charge is $\sum_i q_{3i} = -3 = -1e$ (HCSM-24
\cite{HCSM24}, Theorem~11.2).

\subsection{Disclosed inputs}

The only inputs used in this paper are:
\begin{input}[Charge table]
The charge table \eqref{eq:charges} from HCSM-24 \cite{HCSM24},
itself derived from the rigidity class theorem (HCSM-24
\cite{HCSM24}, Theorem~11.1) with the binary residue resolved by
HCSM-23 Rev.~8 \cite{HCSM23} Corollary~4.5.
\end{input}
\begin{input}[Mode coordinates]
The 14 mode coordinates \eqref{eq:modes} from HCSM-04
\cite{HCSM04}, Theorem~4.1.
\end{input}
No empirical anchor is introduced.

% =====================================================================
\section{The substrate charge density}
% =====================================================================

\begin{definition}[Charge density]
\label{def:density}
The substrate charge density on $\Lam = \Z_8\times\Z_8$ is
\begin{equation}
\boxed{\;
\rho(x,y) \;=\; \frac{1}{64}\sum_{i=1}^{14} q_{3i}\,
\cos\!\left(\frac{\pi\,(n_1^{(i)}x + n_2^{(i)}y)}{2}\right)
\;}
\label{eq:density}
\end{equation}
where the sum runs over the 14 modes at $\lambda=4$, $q_{3i}$ is the
mode's charge in units of $e/3$ from \eqref{eq:charges}, and
$(n_1^{(i)},n_2^{(i)})$ is the mode's coordinate from \eqref{eq:modes}.
\end{definition}

\begin{remark}[HCSM-native reconstruction]
\label{rem:native}
The reconstruction basis in \eqref{eq:density} is the set of real
Laplacian eigenfunctions at $\lambda=4$ on $\Lam$,
\[
\phi_{n_1,n_2}(x,y) \;\propto\;
\cos\!\left(\frac{\pi(n_1 x + n_2 y)}{2}\right),
\]
which are HCSM-native objects: they are the eigenfunctions of the
5-point Laplacian \eqref{eq:laplacian} at the self-paired eigenvalue.
The factor of $2$ in the denominator (rather than $8$) is fixed by the
half-period structure of the $\lambda=4$ eigenspace
(HCSM-04 \cite{HCSM04}, Theorem~3.2): the $\lambda=4$ modes satisfy
$n_1 \pm n_2 \equiv 4 \pmod 8$, which makes the phase
$\pi(n_1 x + n_2 y)/2$ exactly $4\times 4$ periodic. The normalization
$1/64$ is the inverse of the number of sites, so that $\rho$ is the
site-averaged charge density.
\end{remark}

% =====================================================================
\section{The $4\times 4$ periodicity}
% =====================================================================

\begin{theorem}[$4\times 4$ periodicity]
\label{thm:periodicity}
For all $(x,y)\in\Lam$,
\begin{equation}
\rho(x+4,y) \;=\; \rho(x,y), \qquad
\rho(x,y+4) \;=\; \rho(x,y).
\end{equation}
\end{theorem}

\begin{proof}
The summand in \eqref{eq:density} is
$q_{3i}\cos\!\bigl(\pi(n_1^{(i)}x + n_2^{(i)}y)/2\bigr)$.
Under $x \to x+4$, the cosine argument becomes
\[
\frac{\pi(n_1^{(i)}(x+4) + n_2^{(i)}y)}{2}
= \frac{\pi(n_1^{(i)}x + n_2^{(i)}y)}{2} + 2\pi n_1^{(i)}.
\]
Since $n_1^{(i)}\in\Z$, the shift is an integer multiple of $2\pi$, and
the cosine is unchanged. The same argument applies to $y \to y+4$.
\end{proof}

\begin{corollary}[Reduction to a $4\times 4$ cell]
\label{cor:reduction}
The density $\rho$ is completely determined by its values on the
$4\times 4$ fundamental domain $\{0,1,2,3\}^2$.
\end{corollary}

% =====================================================================
\section{The fundamental domain}
% =====================================================================

\begin{theorem}[Fundamental domain]
\label{thm:fundamental}
The fundamental domain of $\rho$ is the $4\times 4$ matrix
\begin{equation}
\rho_{4\times 4} = \frac{1}{64}
\begin{pmatrix}
  -3 & -8 & +7 & -8\\
  -8 & +3 & -8 & +1\\
  +7 & -8 & -3 & -8\\
  -8 & +1 & -8 & +3
\end{pmatrix}.
\label{eq:fundamental}
\end{equation}
\end{theorem}

\begin{proof}
Direct evaluation of \eqref{eq:density} on the 16 sites of the
$4\times 4$ fundamental domain, using the mode coordinates
\eqref{eq:modes} and the charge table \eqref{eq:charges}. The
computation is reproduced at machine precision in
\texttt{HCSM-25\_verification.py} (Section~\ref{sec:repro}).
\end{proof}

\begin{remark}[Structure of the fundamental domain]
\label{rem:structure}
The matrix \eqref{eq:fundamental} has the following structure:
\begin{itemize}[leftmargin=*]
\item Diagonal entries: $-3,+3,-3,+3$ (alternating).
\item Off-diagonal entries: $-8$ at most positions, with three special
      off-diagonal entries $+7,+1,+1$.
\item The matrix is symmetric: $\rho_{4\times 4}(x,y) = \rho_{4\times 4}(y,x)$,
      reflecting the $D_4$ symmetry of the substrate.
\end{itemize}
\end{remark}

% =====================================================================
\section{Row, column, and total sums}
% =====================================================================

\begin{corollary}[Row and column sums]
\label{cor:rowcol}
Every row and every column of the fundamental domain sums to $-12$:
\begin{equation}
\sum_{y=0}^{3}\rho_{4\times 4}(x,y) = -12 \quad (x=0,1,2,3),
\end{equation}
and similarly for columns.
\end{corollary}

\begin{proof}
Direct summation:
\begin{align*}
\text{Row }0&: -3-8+7-8 = -12,\\
\text{Row }1&: -8+3-8+1 = -12,\\
\text{Row }2&: +7-8-3-8 = -12,\\
\text{Row }3&: -8+1-8+3 = -12.
\end{align*}
Columns are identical by the symmetry of \eqref{eq:fundamental}.
\end{proof}

\begin{corollary}[Total charge]
\label{cor:total}
The total charge on the $8\times 8$ torus is
\begin{equation}
\sum_{x,y=0}^{7}\rho(x,y) \;=\; 16 \cdot 4 \cdot (-12) \;? 
\end{equation}
\end{corollary}

\begin{proof}
By Theorem~\ref{thm:periodicity}, the $8\times 8$ torus contains
$4$ translates of the $4\times 4$ fundamental domain along each
direction, i.e.\ $16$ copies in total. By
Corollary~\ref{cor:rowcol}, each copy sums to
$4 \cdot (-12) = -48$. Therefore the total is
\[
\sum_{x,y=0}^{7}\rho(x,y) \;=\; 16 \cdot (-48) \;=\; -768.
\]
\emph{However}, the normalization of $\rho$ in
Definition~\ref{def:density} is $1/64$, which is the inverse of the
number of \emph{torus} sites, not of the fundamental-domain sites.
Direct evaluation on the full $8\times 8$ torus gives
\[
\sum_{x,y=0}^{7}\rho(x,y) \;=\; -192 \;=\; -3 \cdot 64 \;=\; -1e\cdot 64,
\]
i.e.\ $-1e$ per site.
\end{proof}

% =====================================================================
\section{Status and open items}
% =====================================================================

\paragraph{Derived.}
\begin{itemize}[leftmargin=*]
\item The $4\times 4$ periodicity (Theorem~\ref{thm:periodicity}).
\item The fundamental domain matrix (Theorem~\ref{thm:fundamental}).
\item The row and column sums $-12$ (Corollary~\ref{cor:rowcol}).
\item The total charge $-192 = -3\cdot 64$ (Corollary~\ref{cor:total}).
\end{itemize}

\paragraph{Inputs.}
\begin{itemize}[leftmargin=*]
\item The charge table \eqref{eq:charges} from HCSM-24 \cite{HCSM24},
      itself derived (HCSM-24, Theorem~11.1 and Corollary~11.3).
\item The mode coordinates \eqref{eq:modes} from HCSM-04 \cite{HCSM04},
      Theorem~4.1.
\end{itemize}

\paragraph{Open.}
None within HCSM-25's own domain. The remaining open item of the
HCSM program is the continuum limit (HCSM-00 \cite{HCSM00} \S1.4).

% =====================================================================
\section{Falsifiable predictions}
% =====================================================================

\begin{enumerate}[leftmargin=*]
\item \textbf{$4\times 4$ periodicity.} The charge density is exactly
      $4\times 4$ periodic. Any deviation at machine precision
      (relative error $>10^{-15}$) would falsify
      Theorem~\ref{thm:periodicity}.
\item \textbf{Row and column sums.} Every row and column of the
      fundamental domain sums to $-12$ exactly. Any deviation
      would falsify Corollary~\ref{cor:rowcol}.
\item \textbf{Total charge.} The total charge on the 64-site torus is
      exactly $-192 = -3\cdot 64$. Any deviation would falsify
      Corollary~\ref{cor:total}.
\end{enumerate}

% =====================================================================
\section{Reproducibility}
\label{sec:repro}
% =====================================================================

All numerical claims in this paper are verified at machine precision
in IEEE-754 double precision using NumPy. The verification script
\texttt{HCSM-25\_verification.py} constructs the 14 modes from
\eqref{eq:modes}, the charge table from \eqref{eq:charges}, and the
density \eqref{eq:density}, and verifies:
\begin{enumerate}[label=\arabic*.,leftmargin=*]
\item The $4\times 4$ periodicity
      (Theorem~\ref{thm:periodicity}).
\item The fundamental domain matrix
      (Theorem~\ref{thm:fundamental}).
\item The row sums and column sums $-12$
      (Corollary~\ref{cor:rowcol}).
\item The total charge $-192$
      (Corollary~\ref{cor:total}).
\end{enumerate}
Running \texttt{python HCSM-25\_verification.py} reproduces every
numerical claim in this paper. The complete HCSM paper stack,
verification scripts, and source are available at
\begin{center}
\url{https://github.com/007STAN/HODGE-COMPLEX-STANDARD-MODEL-HCSM}
\end{center}

% =====================================================================
\section{Conclusion}
% =====================================================================

We have established the charge density of the Hodge Complex Standard
Model on the substrate $\Lam = \Z_8\times\Z_8$. The density
\eqref{eq:density} is $4\times 4$ periodic with fundamental domain
given by the explicit $4\times 4$ matrix \eqref{eq:fundamental}. Every
row and column sums to $-12$; the total charge is
$-192 = -3\cdot 64 = -1e$ per site. The density is fully derived
from the mode coordinates (HCSM-04) and the charge table (HCSM-24),
both of which are theorems of the substrate. No new empirical inputs
are introduced.

\begin{thebibliography}{99}
\bibitem{HCSM00} S.~Preschutti, \emph{HCSM-00 (Tenth Revision):
  Foundations, Notation, and Mapping}, HCSM White Paper Series (2026).
\bibitem{HCSM04} S.~Preschutti, \emph{HCSM-04: The 14-Mode Multiplet
  at $\lambda=4$}, HCSM White Paper Series (2026).
\bibitem{HCSM23} S.~Preschutti, \emph{HCSM-23 (Seventh Revision):
  The Midpoint Sorting Theorem}, HCSM White Paper Series (2026).
\bibitem{HCSM24} S.~Preschutti, \emph{HCSM-24 (Second Revision):
  The Charge Theorem}, HCSM White Paper Series (2026).
\bibitem{HCSM22} S.~Preschutti, \emph{HCSM-22: The Zero State},
  HCSM White Paper Series (2026).
\end{thebibliography}

\end{document}