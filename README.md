# artifact

C++ implementation of Chadi Saad's algorithm

Pour compiler :
    make std
pour utiliser la librairie standard
    make sparse
pour utiliser sparsepp

Pour lancer un test de rapidité de la librairie standard (remplacer filename par le nom de fichier et X par la taille de k-mère désirée) :

    /usr/bin/time -f "\t%E real,\t%U user,\t%S sys" ./std -f filename -k X > /dev/null

Pour lancer un test de rapidité de la librairie sparsepp (remplacer filename par le nom de fichier et X par la taille de k-mère désirée) :

    /usr/bin/time -f "\t%E real,\t%U user,\t%S sys" ./sparse -f filename -k X > /dev/null
