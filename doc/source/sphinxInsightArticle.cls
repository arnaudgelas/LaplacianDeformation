%
% sphinxInsightArticle.cls for Insight Journal (http://insight-journal.org)
% articles created with sphinx (http://sphinx.pocoo.org/)
%

\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{sphinxInsightArticle}[2013/06/15 Document class (Insight Article)]

% Pass remaining document options to the parent class.
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{\sphinxdocclass}}
\ProcessOptions\relax

\LoadClass{\sphinxdocclass}

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
