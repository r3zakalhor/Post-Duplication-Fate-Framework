\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{manual}[]
% version 1.0 Olivier Pothier 1999/09/15
% version 1.1 Catherine Lamy 2000/04/18

%---------------------------------------
%---- on se base sur la classe book ----
%---- sur papier a4, 2 faces, 1 col ----  
%---------------------------------------
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{book}}
\ProcessOptions
\LoadClass{book}

%----------------------
%---- les packages ----
%----------------------
\usepackage{placeins}
\usepackage{graphics,graphicx,epsfig}
\usepackage[francais,english,ngerman]{babel}
\usepackage[latin1]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{euscript}
\usepackage{battail} 
\usepackage{fancyhdr}
\usepackage{latexsym}
\usepackage{amsmath}
\usepackage{amssymb}
\usepackage{natbib}
\usepackage{array,tabularx}
\usepackage{float}
\usepackage{url}
\usepackage{footnpag}
\usepackage{wrapfig}
%\usepackage[margin=10pt,font={small,sf},labelfont=bf]{caption}
\usepackage[margin=0pt,font={small,sf},labelfont=bf]{caption}
\usepackage[nottoc, notlof, notlot]{tocbibind}

%------------------------------
%---- la taille de la page ----
%------------------------------
%\voffset 0.0cm
%\hoffset 0.0cm
\textheight 23.5cm
\textwidth 16cm
\topmargin 0.0cm
\oddsidemargin 0cm
\evensidemargin 0cm


% Alter some LaTeX defaults for better treatment of figures:
    % See p.105 of "TeX Unbound" for suggested values.
    % See pp. 199-200 of Lamport's "LaTeX" book for details.
    %   General parameters, for ALL pages:
    \renewcommand{\topfraction}{0.9}	% max fraction of floats at top
    \renewcommand{\bottomfraction}{0.8}	% max fraction of floats at bottom
    %   Parameters for TEXT pages (not float pages):
    \setcounter{topnumber}{2}
    \setcounter{bottomnumber}{2}
    \setcounter{totalnumber}{4}     % 2 may work better
    \setcounter{dbltopnumber}{2}    % for 2-column pages
    \renewcommand{\dbltopfraction}{0.9}	% fit big float above 2-col. text
    \renewcommand{\textfraction}{0.07}	% allow minimal text w. figs
    %   Parameters for FLOAT pages (not text pages):
    \renewcommand{\floatpagefraction}{0.7}	% require fuller float pages
	% N.B.: floatpagefraction MUST be less than topfraction !!
    \renewcommand{\dblfloatpagefraction}{0.7}	% require fuller float pages




%----------------------------------------------------------------
%---- pour la sombre histoire que \chapter* ne fait pas appel
%---- à chaptermark, on a toujours \leftmark qui est en 
%---- majuscules (chiant au possible, on remarque au passage
%---- que T. Oetiker, dans ``The Not So Short Introduction to
%---- Latex 2e'' a laissé tombé...). On definit donc un ersazt
%---- de \chapter*, \chapternonum, qui fait appel à \markboth
%---- avec deux arguments identiques (notons au passage que 
%---- \chapter, lui, comme un con, ne passe qu'un argument, le
%---- premier, à markboth au traver de \chaptermark. Il faudrait
%---- aussi réécrire \chapter, demain peut-être :-))
%----------------------------------------------------------------
\newcommand{\chapternonum}[1]{\chapter*{#1\markboth{#1}{#1}}}

\renewcommand{\thechapter}{\Roman{chapter}} 
\renewcommand{\thesection}{\arabic{section}} 
\renewcommand{\thesubsection}{\thesection.\arabic{subsection}} 


%-----------------------------------
%---- on se charge des headings ----
%-----------------------------------
\pagestyle{fancy}
%\bibliographystyle{fancy}
\bibliographystyle{apalike-fr}
%---- on doit redefinir \chaptermark et \sectionmark pour qu'ils
%---- ne mettent pas tout en majuscule. On vire au passage tous les
%---- tests débiles fait dans leur version originale, et on corrige
%---- le fait que \chaptermark ne passait qu'un seul argument à
%---- markboth, ce qui posait problème dans certains cas, rares il est
%---- vrai.
%---- Je ne veux pas afficher \chaptername dans le header, et on
%---- laisse tout dans la casse voulue par l'auteur... L'argument
%---- passé par \chapter à \chaptermark est son argument propre (à
%---- savoir le nom du chapitre. On fait de même pour section.
\renewcommand{\chaptermark}[1]{\markboth {\thechapter.\ #1}{\thechapter.\ #1}}
\renewcommand{\sectionmark}[1]{\markright{\thesection. \ #1}}
\fancyhf{}
%----       le numéro de la page       ----
%---- à gauche sur les pages de gauche ----
%---- à droite sur les pages de droite ----
\fancyhead[LE,RO]{\thepage}
%---- \leftmark, cad le chapitre, à droite des pages de gauche
\fancyhead[RE]{\textsc{\leftmark}}
%---- \rightmark, cad la section, à gauche des pages de droite
\fancyhead[LO]{\textsc{\rightmark}}
%---- les initiales ----%
%\fancyfoot[RE,LO]{\footnotesize C.K.}
%---- la version ----%
%\fancyfoot[CE,CO]{\footnotesize version 1.0}
%---- la date ----%
%\fancyfoot[LE,RO]{\footnotesize \today}
%---- les barres ----%
\renewcommand{\headrulewidth}{1pt}
%\renewcommand{\footrulewidth}{1pt}
\addtolength{\headheight}{2.5pt}
%---- on redefinit le style plain (\chapter le force pour la
%---- première page de tout chapitre
\fancypagestyle{plain}{
  \fancyhead[RE]{}
  \fancyhead[LO]{}
}


\setlength{\parindent}{0pt}

%--------------------------------------------
%---- définitions de nouvelles commandes ----
%--------------------------------------------
\let\urlorig\url
\renewcommand{\url}[1]{%
   \begin{otherlanguage}{english}\urlorig{#1}\end{otherlanguage}%
}

\DeclareGraphicsExtensions{.eps,.png,.pdf}

% Pour que la page paire, ajoutée afin de commencer les chapitres
% sur une page impaire, soit sans header/footer
\newcommand{\clearemptydoublepage}{\newpage{\pagestyle{empty}%
                                   \cleardoublepage}}

\newenvironment{myList}%
{ \begin{list}%
	{$\bullet$}%
	{\setlength{\labelwidth}{30pt}%
	 \setlength{\leftmargin}{30pt}%
	 \setlength{\itemsep}{\parsep}}}%    
{ \end{list} }


\newenvironment{myQuote}[2]%
    {\begin{list}{}{\leftmargin#1\rightmargin#2}\item{}}%
    {\end{list}}


% pour faire des lignes de tableau + épaisses
%%%% debut macro %%%%
\makeatletter
\def\hlinewd#1{%
\noalign{\ifnum0=`}\fi\hrule \@height #1 %
\futurelet\reserved@a\@xhline}
\makeatother
%%%% fin macro %%%%