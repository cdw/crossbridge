% biophysj.cls

\ProvidesClass{biophysj}

\LoadClass[11pt,titlepage]{article}

\usepackage[margin=1in]{geometry}


% Author formatting

% Headers

% Paragraph formatting

% Bibliography

% Tables 
%% Don't change running headers for the list of tables
\renewcommand\listoftables{%
    \section*{\listfigurename}%
    \@starttoc{lof}}
% Table positioning

% Figures
%% Formatting the list of figures by changing all tables of contents
\renewcommand*{\@starttoc}[1]{% 
  \begingroup 
    \makeatletter 
    \parskip\z@ 
    \@input{\jobname.#1}% 
    \if@filesw 
      \expandafter\newwrite\csname tf@#1\endcsname 
      \immediate\openout \csname tf@#1\endcsname \jobname.#1\relax 
    \fi 
    \@nobreakfalse 
  \endgroup 
} 
%% Figures and caption positioning
%% `Change List of Figures' to `Figure Legends'
\renewcommand{\listfigurename}{Figure Legends}
%% Don't change running headers for the list of figures
\renewcommand\listoffigures{%
    \section*{\listfigurename%
        \markboth{FIGURE LEGENDS}{FIGURE LEGENDS}}%
    \@starttoc{lof}}
%% Remove the colon in `Figure X:'
\renewcommand{\fnum@figure}[1]{\figurename~\thefigure}


