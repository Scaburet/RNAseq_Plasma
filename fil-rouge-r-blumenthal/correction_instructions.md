# Instructions de Correction - Fil Rouge RNASeq Partie II

**Date :** mardi 8 octobre 2025  
**Cours :** M1 MEG - UE5 Analyse de données  
**Phase :** Évaluation par les pairs (Phase 2 de l'atelier Moodle)

## Vue d'ensemble

Cette phase d'évaluation par les pairs permet aux étudiants de se corriger mutuellement en utilisant le fichier de correction officiel `GUP_corr_2025.RData`. Les étudiants doivent vérifier que leurs camarades ont correctement fusionné les données de comptage RNASeq.

## Fichiers de correction fournis

- **`GUP_corr_2025.RData`** : Dataframe de correction officiel
- **`FilRouge_RNASeq_R_partII_2025_correction.ipynb`** : Notebook avec solutions complètes
- **`correction_instructions.md`** : Ce document d'instructions

## Instructions pour les étudiants évaluateurs

### 1. Préparation de l'environnement

```r
# Chargement du fichier de correction dans le même répertoire que le notebook à évaluer
load("GUP_corr_2025.RData")

# Vérification que le fichier de correction est bien chargé
ls()  # GUP_corr doit apparaître dans les objets chargés
```

### 2. Chargement du travail à évaluer

```r
# Chargement du fichier de l'étudiant évalué (remplacer par le bon nom)
load("GUP_Nom1_Nom2.RData")

# Vérification que le GUP de l'étudiant est chargé
ls()  # GUP doit apparaître dans les objets chargés
str(GUP)  # Doit afficher un dataframe
```

### 3. Vérifications de structure

```r
# Vérification des dimensions
cat("Dimensions du travail étudiant:", dim(GUP), "\n")
cat("Dimensions attendues:", dim(GUP_corr), "\n")
cat("Dimensions correctes:", identical(dim(GUP), dim(GUP_corr)), "\n")

# Vérification des noms de colonnes
cat("Noms de colonnes corrects:", identical(sort(colnames(GUP)), sort(colnames(GUP_corr))), "\n")

# Vérification du type de données
cat("Type de données correct:", identical(class(GUP), class(GUP_corr)), "\n")
```

### 4. Comparaison avec la correction

```r
# Harmonisation des noms de lignes (peut différer légèrement)
rownames(GUP) <- rownames(GUP_corr)

# Harmonisation des attributs (différences techniques non graves)
attributes(GUP)$row.names <- as.integer(attributes(GUP)$row.names)

# Test d'identité finale
identical_result <- identical(GUP, GUP_corr)
cat("Dataframes identiques:", identical_result, "\n")
```

### 5. Diagnostics en cas de différences

Si les dataframes ne sont pas identiques, effectuer ces vérifications :

```r
# Vérification des colonnes manquantes ou en trop
missing_cols <- setdiff(colnames(GUP_corr), colnames(GUP))
extra_cols <- setdiff(colnames(GUP), colnames(GUP_corr))
cat("Colonnes manquantes:", paste(missing_cols, collapse=", "), "\n")
cat("Colonnes en trop:", paste(extra_cols, collapse=", "), "\n")

# Vérification des gènes manquants ou en trop
missing_genes <- setdiff(rownames(GUP_corr), rownames(GUP))
extra_genes <- setdiff(rownames(GUP), rownames(GUP_corr))
cat("Gènes manquants:", length(missing_genes), "\n")
cat("Gènes en trop:", length(extra_genes), "\n")

# Vérification des valeurs numériques (si structure correcte)
if (identical(dim(GUP), dim(GUP_corr)) && identical(colnames(GUP), colnames(GUP_corr))) {
  numeric_diff <- sum(abs(GUP - GUP_corr), na.rm = TRUE)
  cat("Différence numérique totale:", numeric_diff, "\n")
}
```

## Grille d'évaluation

### Critères d'évaluation (sur 20 points)

| Critère | Points | Description |
|---------|--------|-------------|
| **Structure du dataframe** | 4 pts | Dimensions correctes (25000 × 32) |
| **Nomenclature des colonnes** | 6 pts | Format [génotype]_[sexe]_[réplicat]_[source] respecté |
| **Fusion des données** | 4 pts | Gènes étudiants et auteurs correctement fusionnés |
| **Gestion des identifiants** | 3 pts | Versions de gènes supprimées, mapping SRX/T#### correct |
| **Sauvegarde** | 2 pts | Fichiers .RData et .txt générés |
| **Qualité du code** | 1 pt | Code propre, commenté, sans erreurs |

### Barème détaillé

#### Structure du dataframe (4 points)
- **4 pts** : Dimensions exactes (25000 × 32)
- **3 pts** : Dimensions proches (±5% de gènes, 32 colonnes)
- **2 pts** : Nombre de colonnes correct, nombre de gènes approximatif
- **1 pt** : Structure générale correcte mais dimensions incorrectes
- **0 pt** : Structure incorrecte ou dataframe non créé

#### Nomenclature des colonnes (6 points)
- **6 pts** : Toutes les colonnes suivent le format exact
- **5 pts** : 90% des colonnes correctes
- **4 pts** : 80% des colonnes correctes
- **3 pts** : 70% des colonnes correctes
- **2 pts** : Format général compris mais erreurs systématiques
- **1 pt** : Tentative de renommage mais format incorrect
- **0 pt** : Aucun renommage ou format complètement incorrect

#### Fusion des données (4 points)
- **4 pts** : Fusion parfaite avec gestion des gènes manquants
- **3 pts** : Fusion correcte mais quelques gènes perdus
- **2 pts** : Fusion partielle, données majoritairement présentes
- **1 pt** : Tentative de fusion mais résultat incomplet
- **0 pt** : Pas de fusion ou fusion échouée

#### Gestion des identifiants (3 points)
- **3 pts** : Versions supprimées ET mapping SRX/T#### correct
- **2 pts** : Une des deux tâches correcte
- **1 pt** : Tentative mais erreurs dans l'exécution
- **0 pt** : Identifiants non traités

#### Sauvegarde (2 points)
- **2 pts** : Fichiers .RData et .txt générés et fonctionnels
- **1 pt** : Un seul format sauvegardé
- **0 pt** : Pas de sauvegarde

#### Qualité du code (1 point)
- **1 pt** : Code propre, commenté, exécutable sans erreur
- **0 pt** : Code avec erreurs ou non commenté

## Messages d'erreur courants et solutions

### Erreur : "object 'GUP' not found"
**Cause :** Le fichier .RData de l'étudiant ne contient pas d'objet nommé GUP  
**Solution :** Vérifier que l'étudiant a bien utilisé `save(GUP, file = "...")` et non `save(autre_nom, file = "...")`

### Erreur : "identical() returns FALSE"
**Cause :** Différences dans les données ou la structure  
**Solution :** Utiliser les diagnostics de la section 5 pour identifier les différences

### Erreur : "dimensions differ"
**Cause :** Nombre de gènes ou d'échantillons incorrect  
**Solution :** Vérifier la fusion des dataframes et le nettoyage des identifiants

### Erreur : "column names differ"
**Cause :** Nomenclature incorrecte  
**Solution :** Vérifier le renommage des colonnes selon le format requis

## Commentaires à fournir aux étudiants

### Commentaires positifs (si travail correct)
- "Excellente maîtrise de la manipulation de dataframes en R"
- "Nomenclature parfaitement respectée"
- "Fusion des données réalisée avec succès"
- "Code propre et bien structuré"

### Commentaires constructifs (si erreurs)
- "Attention à la suppression des versions de gènes (ENSMUSG.X → ENSMUSG)"
- "Vérifier le mapping entre identifiants SRX et T#### via les métadonnées"
- "Revoir la nomenclature des colonnes : [génotype]_[sexe]_[réplicat]_[source]"
- "Utiliser merge(..., all = TRUE) pour conserver tous les gènes"

## Contact et support

En cas de questions sur l'évaluation :
- Consulter le notebook de correction complet
- Poser des questions sur le forum Moodle
- Contacter l'équipe pédagogique

---

**Note importante :** Cette évaluation par les pairs fait partie intégrante de l'apprentissage. Soyez bienveillants mais rigoureux dans vos évaluations.
