\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{SNSthesis}

% ifthenelse for if loops
\RequirePackage{ifthen}

% is the thesis for Bachelors, Masters oder PhD? 
\newcommand{\degreetext}{for the degree of Doctor of Philosophy}
\DeclareOption{phd}{\renewcommand{\degreetext}%
{for the degree of \\[0.1in] {\large Doctor of Philosophy}}}
\DeclareOption{master}{\renewcommand{\degreetext}%
{for the degree of \\ Master of Science}}
\DeclareOption{bachelor}{\renewcommand{\degreetext}%
{for the degree of \\ Bachelor of Science}}

% print a copyright page on the back of the title page?
\newboolean{copyrightpage}
\setboolean{copyrightpage}{false}
\DeclareOption{copyrightpage}{\setboolean{copyrightpage}{true}}

% print "Examiner's copy" on the title page?
\newboolean{corrections}
\setboolean{corrections}{false}
\DeclareOption{corrections}{\setboolean{corrections}{true}}

% use small caps for the title text?
\newboolean{titlesmallcaps}
\setboolean{titlesmallcaps}{false}
\DeclareOption{titlesmallcaps}{\setboolean{titlesmallcaps}{true}}

% use fancyschmancy hyperref coloured links, or nice'n black for printing?
\newboolean{foronline}
\setboolean{foronline}{false}
\DeclareOption{foronline}{\setboolean{foronline}{true}}

% finished processing options
\ProcessOptions\relax

% load the book class using a4paper
\LoadClass[12pt,a4paper]{book}

% give the header a bit more room for fancyhdr below
% otherwise LaTeX will spew on each page
\addtolength{\headheight}{2.5pt}

% first set to zero ... 
\setlength{\oddsidemargin}{-1in}
\setlength{\evensidemargin}{-1in}
\setlength{\topmargin}{-1in}       

% adjust these if printer is off by a bit
\setlength{\hoffset}{0mm}
\setlength{\voffset}{0mm}

% from HDR Thesis Preparation Advice 2008
% margins >= 3.5cm on binding edge and >= 1.5cm on opposite
%         >= 2.0cm on top and bottom 

% NB also that the optimal number of characters per line readability is only 60-70, we're over so we'll be a bit more generous on the evensidemargin

\addtolength{\oddsidemargin}{35mm} 
\addtolength{\evensidemargin}{20mm}
\addtolength{\topmargin}{20mm}

% set up some of the spacing
\setlength{\marginparwidth}{40pt}  
\setlength{\marginparsep}{10pt}
\setlength{\headsep}{0.5in}

% A4 dimensions [mm]: 209.903 x 297.039
\setlength{\textwidth}{21 cm}
\setlength{\textheight}{29.7 cm}

% fix up width
\addtolength{\textwidth}{-\oddsidemargin}
\addtolength{\textwidth}{-\evensidemargin}
% now we've added 2inches in setting up margins
\addtolength{\textwidth}{-2in}

% fix up height
\addtolength{\textheight}{-2\topmargin}
\addtolength{\textheight}{-\headheight}
\addtolength{\textheight}{-\headsep}
\addtolength{\textheight}{-\footskip}
% now we've added 2inches in setting up margins
\addtolength{\textheight}{-2in}

\brokenpenalty=10000


% this stops one figure taking up a whole page and lets more text onto the one page when a figure exists
\renewcommand\floatpagefraction{0.8} 
% default = 0.5

% macros for defining the department, various dates, the title and author of the thesis
\def\department#1{\gdef\@department{#1}}
\def\submitdate#1{\gdef\@submitdate{#1}}
\def\copyrightyear#1{\gdef\@copyrightyear{#1}}
% \author, \title in report
\def\@title{}
\def\@author{}
\def\@submitdate{\ifcase\the\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space \number\the\year}
\def\@copyrightyear{\number\the\year}

%% load the required packages
% fancyhdr for nice, fancy headings
\RequirePackage{fancyhdr}
% ccaption for good caption handling
\RequirePackage{ccaption}
% xspace so that spaces after commands are handled correctly
\RequirePackage{xspace}
 % required for nice pictures
\RequirePackage{graphicx} 
% required to use \ifpdf statements, see end of doc
\RequirePackage{ifpdf}

% This constructs the title page
\def\titlep{
    \thispagestyle{empty}
    % \null\vfill
    \begin{center}
        \ifthenelse{\boolean{titlesmallcaps}}
            {\huge\sc\expandafter{\@title}} % put in small caps
            {\Large\uppercase\expandafter{\@title}} % use standard title
    \end{center}
    \vskip2cm
    %\vfill
    \begin{center}
           \textsc{
           PhD program
            \\[0.1in] \large{Nanoscience}
            \\[0.4in] Academic year
            \\[0.1in] \large{2019/2020}}
            \\
	 \vspace{3cm}
        \includegraphics[width=0.72\columnwidth]{logos/itt.jpg}
           \\[0.1in]  \textsc{\large{Scuola Normale Superiore,\\ Istituto Italiano di Tecnologia}}
           % \\[0.1in]  \textsc{\large{Scuola Normale Superiore}}
    \end{center}

        \vfill
        \vspace{1cm}
         \noindent \textsc{Candidate} \hfill \textsc{Thesis Advisor} \par
        \smallskip
        \noindent {\@author} \hfill {Dr. Sergey Artyukhin} \par
        \smallskip
        \noindent \hfill \textsc{Internal Supervisor}\par
        \smallskip
        \noindent \hfill {Prof. Giuseppe La Rocca} \par
        
    \vfill
    \newpage}

% this defines the signature for the signature page where the standard declaration is put
\long\def\signature#1{%
\begin{center}
\begin{minipage}{4in}
\parindent=0pt
Except where acknowledged in the customary manner, the material
presented in this thesis is, to the best of my knowledge, original and
has not been submitted in whole or part for a degree in any university.
\par
\vspace{.5in}
\hbox to 4in{\hfil\shortstack{\vrule width 3in height 0.4pt\\#1}}
\end{minipage}
\end{center}}

% this is the actual signature page
\def\signaturepage{%
	\vskip.5in
	\signature{\@author}}

% now put in the first three or four pages
\def\titlepage{%
    \pagenumbering{roman} % make sure that front matter is numbered roman
    \ifthenelse{\boolean{foronline}}{\hypersetup{urlcolor=black}}{}	
    \titlep
    \ifthenelse{\boolean{copyrightpage}}{
        \newpage
        \vspace*{\fill}
        \begin{center}
            \copyright{} \@author, \@copyrightyear.
        \end{center}
       
        \normalsize
        \newpage
        \signaturepage}{}
	\ifthenelse{\boolean{foronline}}{\hypersetup{urlcolor=blue}}{} % url link in document will be blue
}

% improved version of caption handling
\usepackage{ccaption}
\captionnamefont{\scshape}
\captionstyle{}
\makeatletter
\renewcommand{\fnum@figure}[1]{\quad\small\textsc{\figurename~\thefigure}:}
\renewcommand{\fnum@table}[1]{\quad\small\textsc{\tablename~\thetable}:}
\renewcommand{\@makecaption}[2]{%
\vskip\abovecaptionskip
\sbox\@tempboxa{#1: #2}%
\ifdim \wd\@tempboxa >\hsize
  \def\baselinestretch{1}\@normalsize
  #1: #2\par
  \def\baselinestretch{1.5}\@normalsize
\else
  \global \@minipagefalse
  \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
\fi
\vskip\belowcaptionskip}
\makeatother

% set the pagestyle to look good
\pagestyle{fancy}

% Fancyhdr stuff
% define how headers are marked, for details, see fancyhdr docs
\renewcommand{\chaptermark}[1]{\markboth{#1}{}}
\renewcommand{\sectionmark}[1]{\markright{\thesection\ #1}}

% (alternatively some stuff done by Martin)
%\renewcommand{\sectionmark}[1]{\markright{\thesection\ #1}}
%\makeatletter
%\renewcommand{\chaptermark}[1]{%
%  \markboth{\if@mainmatter\textit{\chaptername\ \thechapter\  }\fi #1}{}%
%}
%\makeatother


% define where sections, chapters and pagenumbers are put
% see fancyhdr docs for details
% the \nouppercase stops book.cls making the contents, bibliography and index headers from being all in uppercase.
% The options used here are essentially that in Lamport's book, but with small caps for the headings.
\fancyhf{}
\fancyhead[LE,RO]{\nouppercase{\thepage}}
\fancyhead[LO]{\sc \nouppercase{\rightmark}}
\fancyhead[RE]{\sc \nouppercase{\leftmark}}

% (alternatively some stuff done by Martin)
%\fancyhf{}
%\fancyhead[LE,RO]{\nouppercase{\thepage}}
%\fancyhead[LO]{\sc \nouppercase{\rightmark}}
%\fancyhead[RE]{\rm \nouppercase{\leftmark}}


%% other settings required for a thesis
% use 1.5 line spacing for the examiners copy
% single line spacing for the final copy
\ifthenelse{\boolean{corrections}}{%
\renewcommand{\baselinestretch}{1.5}}%
{\renewcommand{\baselinestretch}{1}}

% spell things correctly
\newenvironment{centre}{\begin{center}}{\end{center}}
\newenvironment{itemise}{\begin{itemize}}{\end{itemize}}

%\usepackage{play}
\usepackage[grey,times]{quotchap} % this makes the chapter title look nice and you can insert a quote
%\renewcommand{\chapterheadstartvskip}{\vspace*{1.1\baselineskip}} % changes vspace for quotes
\usepackage{makeidx}

%% set up the bibliography style
\bibliographystyle{SNSthesis}  % SNSthesis bibliography style file, made with makebst

%% optional packages
\usepackage[square,comma,numbers,sort&compress]{natbib}
% this is the natural sciences bibliography citation style package.
% The options here give citations in the text as numbers in square brackets, separated by commas, citations sorted and consecutive citations compressed
% output example: [1,4,12-15]

\usepackage[nottoc]{tocbibind}  
% allows the table of contents, bibliography and index to be added to the table of contents if desired, the option used here specifies that the table of contents is not to be added.
% tocbibind needs to be after natbib otherwise bits of it get trampled.


\ifpdf
%-for more info checkout: http://en.wikibooks.org/wiki/LaTeX/Hyperlinks where you can also find the manual with more options.
	\ifthenelse{\boolean{foronline}}{
          \usepackage[pdftex, plainpages = false, pdfpagelabels, 
                 pdfpagelayout = OneColumn,
                 bookmarks,
                 bookmarksopen = true,
                 bookmarksnumbered = true,
                 breaklinks = true,
                 linktocpage,
                 hyperindex = true,
                 hyperfigures,
		pdfdisplaydoctitle = true,
                pagebackref,
                colorlinks=true,
                linkcolor=blue,
                urlcolor=blue,
                citecolor=red,
                anchorcolor=green
                 ]{hyperref}}{\usepackage[pdftex]{hyperref}}
    	\pdfcompresslevel=9
\else
      \ifthenelse{\boolean{foronline}}{
        \usepackage[ dvips, 
                 bookmarks,
                 bookmarksopen = true,
                 bookmarksnumbered = true,
                 breaklinks = true,
                 linktocpage,
                 hyperindex = false,
                 hyperfigures,
                 pagebackref,
                 colorlinks=true,
                 linkcolor=blue,
                 urlcolor=blue,
                 citecolor=red,
                 anchorcolor=green
                 ]{hyperref}}{\usepackage[dvips]{hyperref}}	
\fi

% put in an index?
%\makeindex
