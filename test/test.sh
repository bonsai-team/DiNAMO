#!/bin/bash

TIMECMD=/usr/bin/time
FILE1=../sequences/small.fasta
FILE2=../sequences/medium.fasta
FILE3=../sequences/chr22.fa
SCRIPTDIR=`dirname "$(readlink -f "$0")"`

if [[ $1 = "exec_time" ]]
then
    yaxis="Secondes d'éxécution"
    graphTitle="Temps d'exécution en fonction de la taille des fichiers et de la taille des k-mères"
    TIMEOPT="-f %e"
else
    yaxis="Mémoire maximum utilisée (en Ko)"
    graphTitle="Mémoire maximum utilisée en fonction de la taille des fichiers et de la taille des k-mères"
    TIMEOPT="-f %M"
fi

#OTZ
header="<!DOCTYPE html>\n<html>\n<head>\n\t<meta http-equiv=\"content-type\" content=\"text/html; charset=UTF-8\">\n\n\t<title></title>\n</head>\n<body>\n\t<script src=\"$SCRIPTDIR/lib/highcharts.js\"></script>\n\t<script src=\"$SCRIPTDIR/lib/exporting.js\"></script>\n\n\t<div id=\"container\" style=\"min-width: 400px; height: 400px; margin: 0 auto\"></div>\n\n\t<script type='text/javascript'>\n\nHighcharts.chart('container', {\n\tchart: {\n\t\ttype: 'column'\n\t},\n\ttitle: {\n\t\ttext: \"$graphTitle\"\n\t},\n\txAxis: {\n\t\tcategories: [\n\t\t\t'Petit fichier',\n\t\t\t'Fichier moyen',\n\t\t\t'Grand fichier'\n\t\t]\n\t},\n\tyAxis: [{\n\t\tmin: 0,\n\t\ttitle: {\n\t\t\ttext: \"$yaxis\"\n\t\t}\n\t}],\n\tplotOptions: {\n\t\tseries: {\n\t\t\tgrouping: false,\n\t\t\tshadow: false\n\t\t}\n\t},\n\tseries: [{\n\t\tname: 'Librairie Standard',\n\t\tid: 'S'\n\t}, {\n\t\tname: 'Sparsepp',\n\t\tid: 'C'\n\t}"

footer="]\n});\n\n</script>\n\n</body>\n\n</html>"

# réinitialisation de la page HTML
rm -f graph.html && echo -en $header > graph.html

make sparse
make std

for (( i=4; i<=12; i+=2))
do
    std_title="Standard/KM=$i"
    std_color="#E74C3C"
    std_id="S"
    pointPlacement=`echo "-0.40 + ($i - 4)*0.1" | bc`
    std_pointWidth="26"
    echo -en ",{\n\t\tname: \"$std_title\",\n\t\tcolor: \"$std_color\",\n\t\tdata: [" >> graph.html
    echo -en "$($TIMECMD "$TIMEOPT" ./std -f $FILE1 -k $i 2>&1 >/dev/null), " >> graph.html
    echo -en "$($TIMECMD "$TIMEOPT" ./std -f $FILE2 -k $i 2>&1 >/dev/null), " >> graph.html
    echo -en "$($TIMECMD "$TIMEOPT" ./std -f $FILE3 -k $i 2>&1 >/dev/null)],\n" >> graph.html
    echo -en "\t\tpointPlacement: $pointPlacement,\n\t\tpointWidth: $std_pointWidth,\n\t\tlinkedTo: \"$std_id\"\n\t}" >> graph.html

    sparse_title="Sparsepp/KM=$i"
    sparse_color="#3498DB"
    sparse_id="C"
    sparse_pointWidth="16"
    echo -en ",{\n\t\tname: \"$sparse_title\",\n\t\tcolor: \"$sparse_color\",\n\t\tdata: [" >> graph.html
    echo -en "$($TIMECMD "$TIMEOPT" ./sparse -f $FILE1 -k $i 2>&1 >/dev/null), " >> graph.html
    echo -en "$($TIMECMD "$TIMEOPT" ./sparse -f $FILE2 -k $i 2>&1 >/dev/null), " >> graph.html
    echo -en "$($TIMECMD "$TIMEOPT" ./sparse -f $FILE3 -k $i 2>&1 >/dev/null)],\n" >> graph.html
    echo -en "\t\tpointPlacement: $pointPlacement,\n\t\tpointWidth: $sparse_pointWidth,\n\t\tlinkedTo: \"$sparse_id\"\n\t}" >> graph.html

done

echo -e $footer >> graph.html

make clean
