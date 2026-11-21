# Build from repository root: latexmk papers/nc-p3.tex
# research/plan_ja.tex and research/notes_ja.tex use the same build directory.
@default_files = ('papers/*.tex');
$do_cd = 1;
$pdf_mode = 4;
$lualatex = 'lualatex -synctex=1 -interaction=nonstopmode -halt-on-error -file-line-error %O %S';
$out_dir = '../build';
$aux_dir = $out_dir;
$max_repeat = 6;
