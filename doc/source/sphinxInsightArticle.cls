%
% sphinxInsightArticle.cls for Insight Journal (http://insight-journal.org)
% articles created with sphinx (http://sphinx.pocoo.org/)
%

\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{sphinxInsightArticle}[2013/06/15 Document class (Insight Article)]

%
%  Change this to say a4paper instead of letterpaper if you want A4.  These
%  are the latex defaults.
%
\newcommand{\itk@paper}{letterpaper}
\newcommand{\itk@ptsize}{11pt}

%
%  Commands for image and figure captions
%
\newcommand{\itkcaption}[2][]{\caption[{#1}]{\small \textsf{{#2}} \normalsize}}
\newcommand{\itkpiccaption}[2][]{\piccaption[{#1}]{\small \textsf{{#2}} \normalsize}}

%  These set up the fonts for the documents.
%
%  The "times" package makes the default font the PostScript Times
%  font, which makes for smaller PostScript and a font that more people
%  like.
%
\RequirePackage{times}\typeout{Using Times instead of Computer Modern.}

% Change the options here to get a different set of basic options,  This
% is where to add things like "a4paper" or "10pt".
%
\LoadClass[\itk@paper,\itk@ptsize]{article}

\setcounter{secnumdepth}{1}

% Optional packages:
%
% If processing of these documents fails at your TeX installation,
% these may be commented out (independently) to make things work.
% These are both supplied with the current version of the teTeX
% distribution.
%
% The "fancyhdr" package makes nicer page footers reasonable to
% implement, and is used to put the chapter and section information in
% the footers.
%
\RequirePackage{fancyhdr}\typeout{Using fancier footers than usual.}

% Adding the package ifpdf in order to deal with the cases where
% pdflatex is used instead of latex.
\RequirePackage{ifpdf}\typeout{Using ifpdf for the template}

% Required package:
%
% This gives us all the Insight-specific markup that we really want.
% This should come last.  Do not change this.
%
\RequirePackage{sphinxInsightJournal}

%--------------------------------------------- from sphinx below
% This is needed to get the width of the section # area wide enough in the
% library reference.  Doing it here keeps it the same for all the manuals.
%
\renewcommand*\l@section{\@dottedtocline{1}{1.5em}{2.6em}}
\renewcommand*\l@subsection{\@dottedtocline{2}{4.1em}{3.5em}}

% Fix the bibliography environment to add an entry to the Table of
% Contents.
% For a report document class this environment is a chapter.
\let\py@OldThebibliography=\thebibliography
\renewcommand{\thebibliography}[1]{
  \cleardoublepage
  \phantomsection
  \py@OldThebibliography{1}
  \addcontentsline{toc}{chapter}{\bibname}
}
