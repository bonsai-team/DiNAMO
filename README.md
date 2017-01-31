# artifact

C++ implementation of Chadi Saad's algorithm

Pour compiler :
make

Pour lancer un test de rapidité (remplacer filename par le nom de fichier et X par la taille de k-mère désirée) :

    /usr/bin/time -f "\t%E real,\t%U user,\t%S sys" ./hash -f filename -k X > /dev/null

Pour le moment je n'ai pas réussi à faire un switch interactif entre la lib standard et la lib optimisée, mais on peut le faire facilement en bougeant le commentaire l.85/86 :

    sparse_hash_map<string, int> encounters;
    //map<string, int> encounters;

