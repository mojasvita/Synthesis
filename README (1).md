# Synteny / protein homology pipeline

Минимальный анализ одного гена *D. melanogaster* (`FBgn...`) в сравнительных геномах из HAL.

## Запуск

Из папки `hack` в подготовленном окружении (`phylo-ursidae`):

```bash
conda activate phylo-ursidae

GENE=FBgn0035694
python 01_prepare_query.py "$GENE"
python 02_find_syntenic_loci.py "$GENE"
python 03_blast_and_report.py "$GENE" --evalue 1e-5 --min-qcov 70
python 04_spaln_uncertain_loci.py "$GENE" --flank 10000
python 05_species_tree.py "$GENE"
```

Нужны `Acalyptratae20.hal`, файлы GFF3 / Protein FASTA в `reference_data/`, дерево `reference_data/phylogeny/species_tree.nwk`, а также установленные `hal2fasta`, `halLiftover`, `blastp`, `spaln`, Python с Biopython и Matplotlib. Пути и параметры можно уточнить через `python <script>.py --help`.

## Результаты

В `results/<GENE>/`:

- `locus_states.tsv` — синтенические локусы, пересечения с аннотированными генами, BLASTP-hit и статусы.
- `spaln/spaln_run_summary.tsv` и `spaln/*/spaln_alignment.txt` — результаты Spaln **только для локусов, выбранных его фильтром** (если таких локусов нет, файлов может не быть).
- `locus_states_with_spaln.tsv` — объединённые статусы (без перезаписи исходного TSV).
- `species_status.tsv` — один статус на вид.
- `species_tree.png`, `species_tree.pdf`, `species_tree.svg` — визуализация на видовом дереве.

## Ограничения

- Параметры BLASTP выше (`E-value < 1e-5`, query coverage ≥ 70%) — **рабочие пороги**, не подтверждение ортологии и не обязательно пороги статьи. Отдельного фильтра минимальной длины 40 aa в скрипте **нет**.
- BLASTP сравнивает белок с белками **синтенического локуса**, а не со всеми белками протеома; E-value зависит от размера базы поиска.
- Spaln проверяет геномный участок с флангами ±10 kb, но найденное выравнивание **не считается автоматически достоверным**. Скрипт 04 обрабатывает `annotated_without_supported_homology` и `annotated_without_blast_hit`, но **не** `no_coding_annotation_needs_genewise`.
- `uncertain` означает **недостаток доказательств**, а не отсутствие гена или noncoding-состояние. `not_assessed` — вид не обработан.
- Species tree отображает статусы на готовой топологии HAL; само по себе дерево **не доказывает de novo происхождение**. Для этого нужны дополнительные проверки предковых некодирующих локусов, отдалённых гомологов и паралогов.
