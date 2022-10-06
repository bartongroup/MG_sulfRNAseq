
rsync -rvm --include='*/' --include='.htaccess' --include='doc/report.html' --include='tab/*' --exclude='*' . cluster:/cluster/gjb_lab/mgierlinski/public_html/covid_sulf_rnaseq