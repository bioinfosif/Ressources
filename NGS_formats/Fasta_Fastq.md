# Les formats de fichiers NGS : FASTA et FASTQ

## Introduction

Le séquençage à haut débit (NGS - Next Generation Sequencing) génère d'énormes quantités de données qui doivent être stockées et analysées efficacement. Deux formats de fichiers sont particulièrement importants dans ce domaine : FASTA et FASTQ. Ces formats servent à stocker des séquences biologiques (ADN, ARN ou protéines) et, dans le cas du FASTQ, des informations sur la qualité de chaque base séquencée.

## Le format FASTA

### Structure du format FASTA

Le format FASTA est un format texte simple utilisé pour représenter des séquences de nucléotides ou d'acides aminés. Chaque séquence dans un fichier FASTA comprend deux parties principales :

1. **Une ligne d'en-tête** : Elle commence toujours par le symbole ">" suivi d'un identifiant et optionnellement d'une description.
2. **La séquence elle-même** : Représentée par une série de lettres correspondant aux nucléotides (A, T, G, C pour l'ADN) ou aux acides aminés.

### Exemple de fichier FASTA

```
>Seq1 description de la première séquence
ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
>Seq2 description de la deuxième séquence
CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC
GATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
```

### Caractéristiques du format FASTA

- **Format simple et lisible** par l'humain
- **Économe en espace** par rapport à d'autres formats
- **Largement utilisé** pour stocker tout type de séquences biologiques
- **Ne contient pas d'information sur la qualité** des séquences

### Applications courantes

- Stockage de séquences de référence (génomes, protéomes)
- Entrée/sortie pour de nombreux outils bioinformatiques
- Partage de séquences entre chercheurs
- Base de données génomiques


## Le format FASTQ

### Structure du format FASTQ

Le format FASTQ étend le format FASTA en ajoutant des informations sur la qualité du séquençage pour chaque base. Un fichier FASTQ contient quatre lignes par séquence :

1. **Ligne d'identifiant** : Commence par "@" suivi d'un identifiant unique et éventuellement d'une description
2. **Ligne de séquence** : Les nucléotides bruts (A, T, G, C, N)
3. **Ligne séparateur** : Généralement un simple "+"
4. **Ligne de qualité** : Encodage ASCII des scores de qualité pour chaque base de la séquence

### Exemple de fichier FASTQ

```
@SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=36
GGGTGATGGCCGCTGCCGATGGCGTCAAATCCCACC
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IG9IC
@SRR001666.2 071112_SLXA-EAS1_s_7:5:1:801:338 length=36
CAGGCGCCCGACCGCGTCTTTGACTTTCTTGAAAAA
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII6IBI
```

### Scores de qualité Phred

Les scores de qualité représentent la probabilité d'erreur pour chaque base séquencée. Ils sont encodés en caractères ASCII pour économiser de l'espace.

| Score Phred | Probabilité d'erreur | Précision |
|-------------|----------------------|-----------|
| 10          | 1 sur 10             | 90%       |
| 20          | 1 sur 100            | 99%       |
| 30          | 1 sur 1,000          | 99.9%     |
| 40          | 1 sur 10,000         | 99.99%    |

La formule du score Phred : Q = -10 × log₁₀(P)
où P est la probabilité d'erreur.

### Encodages des scores de qualité

Plusieurs encodages existent pour les scores de qualité :

- **Phred+33** (Sanger) : ASCII 33-126
- **Phred+64** (Illumina 1.3+) : ASCII 64-126
- **Phred+33 v2** (Illumina 1.8+) : ASCII 33-74


### Caractéristiques du format FASTQ

- **Contient les données brutes** de séquençage
- **Inclut les scores de qualité** pour chaque base
- **Format standard** pour les données de sortie des séquenceurs NGS
- **Taille plus importante** que le format FASTA en raison des informations de qualité

## Comparaison FASTA vs FASTQ

| Caractéristique | FASTA | FASTQ |
|-----------------|-------|-------|
| Ligne d'en-tête | Commence par ">" | Commence par "@" |
| Information de qualité | Non | Oui |
| Nombre de lignes par séquence | 2 | 4 |
| Taille relative | Plus petite | Plus grande |
| Utilisation principale | Séquences de référence | Données brutes de séquençage |

## Traitement des fichiers FASTA et FASTQ

### Outils courants

1. **FastQC** : Contrôle qualité des données FASTQ
2. **Trimmomatic** : Filtrage et nettoyage des séquences FASTQ
3. **BLAST** : Recherche de similarités dans les séquences FASTA
4. **Bowtie/BWA** : Alignement de séquences FASTQ sur un génome de référence FASTA
5. **SAMtools** : Manipulation des fichiers d'alignement générés à partir de FASTQ

### Flux de travail typique en NGS

1. Séquençage → fichiers FASTQ
2. Contrôle qualité et filtrage des fichiers FASTQ
3. Alignement sur un génome de référence (FASTA)
4. Analyse en aval (variant calling, RNA-Seq, etc.)


## Exemples pratiques

### Conversion FASTQ → FASTA

La conversion d'un fichier FASTQ en FASTA implique de supprimer les informations de qualité et de changer le caractère d'en-tête de "@" à ">".

**Script Python simple pour la conversion :**

```python
def fastq_to_fasta(fastq_file, fasta_file):
    with open(fastq_file, 'r') as in_file, open(fasta_file, 'w') as out_file:
        line_count = 0
        for line in in_file:
            line_count += 1
            if line_count % 4 == 1:  # Ligne d'identifiant
                out_file.write('>' + line[1:])  # Remplace @ par >
            elif line_count % 4 == 2:  # Ligne de séquence
                out_file.write(line)
```

### Extraction des scores de qualité moyens d'un fichier FASTQ

```python
def average_quality_score(fastq_file):
    with open(fastq_file, 'r') as in_file:
        line_count = 0
        total_quality = 0
        total_bases = 0
        
        for line in in_file:
            line_count += 1
            if line_count % 4 == 0:  # Ligne de qualité
                for char in line.strip():
                    # Convertit le caractère ASCII en score Phred (pour Phred+33)
                    quality = ord(char) - 33
                    total_quality += quality
                    total_bases += 1
                    
        return total_quality / total_bases if total_bases > 0 else 0
```

## Conclusion

Les formats FASTA et FASTQ sont fondamentaux dans l'analyse des données de séquençage moderne. Bien que simples dans leur conception, ils permettent de stocker efficacement d'énormes quantités d'informations génomiques.

- **FASTA** est idéal pour les séquences de référence et les analyses qui ne nécessitent pas d'information sur la qualité des bases.
- **FASTQ** est essentiel pour les données brutes de séquençage où la qualité de chaque base doit être prise en compte pour les analyses en aval.

La maîtrise de ces formats et des outils associés est indispensable pour tout bioinformaticien travaillant dans le domaine du NGS.

## Ressources supplémentaires

- [Documentation de NCBI sur le format FASTA](https://ncbi.nlm.nih.gov/blast/fasta.shtml)
- [Spécification du format FASTQ](https://en.wikipedia.org/wiki/FASTQ_format)
- [Galaxy Project - Tutoriels sur l'analyse NGS](https://galaxyproject.org/learn/)
- [Bioconductor - Packages R pour l'analyse de données NGS](https://bioconductor.org/)
