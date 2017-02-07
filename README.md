# artifact

C++ implementation of Chadi Saad's algorithm

Pour compiler :

    make std

pour utiliser la librairie standard

    make sparse

pour utiliser sparsepp

Pour lancer un test de rapidité de la librairie standard (remplacer filename par le nom de fichier et X par la taille de k-mère désirée) :

    /usr/bin/time -f "\tElapsed time : %Es \t\t Maximum memory used : %MkB" ./std -f filename -k X > /dev/null

Pour lancer un test de rapidité de la librairie sparsepp (remplacer filename par le nom de fichier et X par la taille de k-mère désirée) :

    /usr/bin/time -f "\tElapsed time : %Es \t\t Maximum memory used : %MkB" ./sparse -f filename -k X > /dev/null
