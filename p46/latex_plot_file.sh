
old_dir=$PWD

cd $plot_dir

latex_plot_file=plot_$arg

cat << EOF > $latex_plot_file.tex 
\documentclass{report}
\usepackage{graphicx}
\begin{document}
\begin{figure}[ht]
\includegraphics[width=$latex_plot_width]{$arg.eps}
\end{figure}
\caption{$latex_plot_caption}
\end{figure}
\end{document}
EOF

latex $latex_plot_file 
dvitopdf $latex_plot_file

cd $old_dir
