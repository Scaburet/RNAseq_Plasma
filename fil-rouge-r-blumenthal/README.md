# Fil Rouge RNASeq - Partie II : Manipulation de Dataframes en R

**Date :** mardi 8 octobre 2025  
**Cours :** M1 MEG - UE5 Analyse de données  
**Enseignants :** Claire Vandiedonck, Sandrine Caburet

## Vue d'ensemble

Ce module pédagogique porte sur la manipulation de dataframes en R dans le contexte de l'analyse de données RNASeq. Les étudiants apprennent à fusionner des données de comptage provenant de deux pipelines bioinformatiques différents.

## Contexte scientifique

- **Organisme :** Souris (*Mus musculus*)
- **Région d'intérêt :** Chromosome 7 (région impliquée dans les troubles du spectre autistique)
- **Conditions :** Délétion, duplication, wild type
- **Données :** Bulk RNASeq, 16 échantillons
- **Référence :** Blumenthal et al., 2014

## Structure du module

### Fichiers pour les étudiants

1. **`FilRouge_RNASeq_R_partII_2025.ipynb`**
   - Notebook principal avec exercices guidés
   - 22 cellules (markdown + code R)
   - Instructions détaillées et exemples progressifs

### Fichiers de correction

2. **`FilRouge_RNASeq_R_partII_2025_correction.ipynb`**
   - Solutions complètes avec explications détaillées
   - 21 cellules avec vérifications et diagnostics
   - Gestion d'erreurs et cas particuliers

3. **`GUP_corr_2025.RData`**
   - Dataframe de correction officiel pour l'évaluation par les pairs
   - 25,000 gènes × 32 échantillons
   - Nomenclature correcte : `[génotype]_[sexe]_[réplicat]_[source]`

4. **`correction_instructions.md`**
   - Instructions détaillées pour l'évaluation par les pairs
   - Grille d'évaluation sur 20 points
   - Diagnostics et résolution d'erreurs

## Objectifs pédagogiques

### Compétences techniques
- Lecture de fichiers de données avec différents formats
- Manipulation de dataframes : renommage, fusion, extraction
- Gestion des identifiants et versions de gènes
- Sauvegarde dans différents formats (.RData, .txt)

### Compétences transversales
- Résolution de problèmes d'intégration de données
- Compréhension des différences entre pipelines bioinformatiques
- Travail collaboratif via l'évaluation par les pairs

## Défis techniques abordés

1. **Identifiants hétérogènes**
   - Étudiants : SRX1589825, SRX1589826, etc.
   - Auteurs : T1979_ctx, T1980_ctx, etc.
   - Solution : Mapping via fichier de métadonnées

2. **Versions de gènes**
   - Étudiants : ENSMUSG00000051951.6 (avec version)
   - Auteurs : ENSMUSG00000051951 (sans version)
   - Solution : Suppression des versions avec regex

3. **Nomenclature complexe**
   - Format requis : `[génotype]_[sexe]_[réplicat]_[source]`
   - Exemple : `wtdel_f_1_s` (wild type deletion, female, replicate 1, students)

4. **Fusion de données**
   - Gènes présents dans un seul jeu de données
   - Solution : Outer join avec remplacement NA → 0

## Données d'entrée

### Fichiers sur le serveur Adenine
- **Comptages étudiants :** `/srv/data/meg-m1-a2b/blumenthal-2014/wholeData/4-counts/16samples_counts_matrix.txt`
- **Comptages auteurs :** `/srv/data/meg-m1-a2b/blumenthal-2014/wholeData/4-counts/GSE76872_cortexGeneCounts_authors.csv`
- **Métadonnées :** `/srv/data/meg-m1-a2b/blumenthal-2014/wholeData/4-counts/SraRunTable-Blumenthal2014-cortex.tsv`

### Structure attendue du résultat final
```
GUP dataframe : 25,000 gènes × 32 colonnes
├── 16 colonnes étudiants : *_s (students)
├── 16 colonnes auteurs : *_a (authors)
└── Nomenclature : wtdel_f_1_s, wtdel_f_1_a, wtdel_f_2_s, ...
```

## Environnement technique

- **Serveur :** adenine.plasmabio.org
- **Ressources :** 5 CPUs, 10 GB RAM par étudiant
- **R version :** 4.x avec tidyverse
- **Format :** Jupyter notebook avec kernel R

## Évaluation

### Modalités
- **Type :** Travail à la maison évalué
- **Méthode :** Évaluation par les pairs (Phase 2 atelier Moodle)
- **Barème :** 20 points
- **Livrables :** Notebook complété + fichiers GUP_Nom1_Nom2.RData et .txt

### Critères principaux
1. Structure du dataframe (4 pts)
2. Nomenclature des colonnes (6 pts)
3. Fusion des données (4 pts)
4. Gestion des identifiants (3 pts)
5. Sauvegarde (2 pts)
6. Qualité du code (1 pt)

## Installation et utilisation

### Pour les enseignants
1. Déployer les notebooks sur le serveur Adenine
2. Vérifier l'accès aux fichiers de données
3. Distribuer les instructions d'évaluation

### Pour les étudiants
1. Ouvrir `FilRouge_RNASeq_R_partII_2025.ipynb` sur Adenine
2. Suivre les instructions étape par étape
3. Sauvegarder les résultats avec vos noms de famille
4. Participer à l'évaluation par les pairs

## Adaptation depuis l'année précédente

### Changements principaux
- **Données :** CMH humain → Chromosome 7 souris
- **Échantillons :** 6 → 16
- **Conditions :** 3 lignées cellulaires → 3 génotypes × 2 sexes
- **Pipeline :** Génomes personnalisés → Génome de référence unique

### Éléments conservés
- Structure pédagogique progressive
- Évaluation par les pairs
- Niveau de difficulté adapté au M1
- Format Jupyter notebook

## Support et contact

- **Questions techniques :** Forum Moodle
- **Problèmes d'accès :** Équipe technique Adenine
- **Contenu pédagogique :** Claire Vandiedonck, Sandrine Caburet

---

*Ce module s'inscrit dans le cadre de la nouvelle maquette 2025 du M1 MEG et contribue à l'évaluation continue de l'UE5.*
