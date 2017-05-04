# DiNAMO

[![Build Status](https://travis-ci.org/bonsai-team/DiNAMO.svg?branch=master)](https://travis-ci.org/bonsai-team/DiNAMO)

C++ implementation of Chadi Saad's algorithm

## Compilation

Pour compiler :

    make

## Exécution

Pour exécuter le programme vous devez fournir les options suivantes :

  * -nf *chemin/vers/votre/fichier/negatif*
  * -pf *chemin/vers/votre/fichier/positif*
  * -l  *taille\_de\_motif\_désirée*

  Une ligne de commande typique pourrait être :

        ./artifact -nf sequences/controle.fa -pf sequences/signal.fa -l 6

### Options

Les options suivantes sont aussi disponibles :

  * limiter le nombre de positions qui pourront être dégénérées

        -d limite

  * ne compter que les motifs qui sont à une certaine position dans la séquence

        -p position_du_motif

  (Veuillez noter que le comptage se fait depuis la fin de la séquence, ainsi la position 0 correspond au dernier motif de chaque séquence)
