
# artifact

C++ implementation of Chadi Saad's algorithm

## Compilation

Pour compiler :

    make

## Exécution

Pour exécuter le programme vous devez fournir les options suivantes :

  * -f *chemin/vers/votre/fichier*
  * -k *taille\_de\_motif\_désirée*

  Une ligne de commande typique pourrait être :

        ./artifact -f sequences/chr22.fa -k 6

Les options suivantes sont aussi disponibles :

  * limiter le nombre de positions qui pourront être dégénérées

        -d limite

  * ne compter que les motifs qui sont à une certaine position dans la séquence

  (Veuillez noter que le comptage se fait depuis la fin de la séquence, ainsi la position 0 correspond au dernier motif)

        -p position_du_motif


# Tests

Rendez-vous dans le répertoire test

## Test de rapidité d'éxecution

Pour lancer un test de rapidité entre la librairie standard et sparse++ (résultat produit dans graph.html) :

    make exectime

## Test d'empreinte mémoire

Pour lancer un test de mémoire entre librairie standard et sparse++ (résultat produit dans graph.html) :

    make maxmemory
