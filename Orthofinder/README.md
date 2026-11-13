# Вариант работы с Orthofinder

**Исходные данные**

Для сравнительного анализа использовались белковые протеомы 28 видов.

В набор входили:

17 видов рода Drosophila;

другие двукрылые насекомые;

представители других отрядов насекомых;

внешний сравнительный таксон Daphnia magna.

В частности, в анализ включены:

Drosophila:
D. melanogaster, D. simulans, D. sechellia, D. yakuba, D. erecta, D. ananassae, D. biarmipes, D. elegans, D. eugracilis, D. ficusphila, D. kikkawai, D. miranda, D. persimilis, D. pseudoobscura, D. rhopaloa, D. takahashii, D. willistoni.

Другие таксоны:
Aedes aegypti, Anopheles gambiae, Apis florea, Apis mellifera, Bactrocera dorsalis, Bactrocera latifrons, Bombyx mori, Daphnia magna, Monomorium pharaonis, Musca domestica, Scaptodrosophila lebanonensis.

Общий набор содержал 398 017 белковых последовательностей.

**3. Использованные программы**

OrthoFinder 3.2

Поиск ортогрупп, ортологов, событий дупликации и построение эволюционных деревьев генов

DIAMOND 2.2.8.188

Поиск сходства белковых последовательностей

MAFFT

Множественное выравнивание белковых последовательностей

FastTree

Построение генетических деревьев

Python / Biopython

Обработка FASTA, таблиц OrthoFinder и результатов DIAMOND

**Запуск Orthofinder**

`orthofinder \

    -f drosophila_proteomes/proteomes_1isoform \

    -S diamond_ultra_sens \

    -M msa \

    -A mafft \

    -T fasttree \

    -t 16 \

    -a 8
`
OrthoFinder использовал DIAMOND для поиска сходства последовательностей, MAFFT для выравнивания и FastTree для построения деревьев.

## Основные результаты OrthoFinder

Для последующей классификации использовались следующие файлы:

| Файл | Использование |
|------|---------------|
| `Orthogroups/Orthogroups.tsv` | Состав ортогрупп и принадлежность белков к ним |
| `Orthogroups/Orthogroups_UnassignedGenes.tsv` | Белки, не распределённые по ортогруппам |
| `Orthogroups/Orthogroups.GeneCount.tsv` | Количество белков каждого вида в ортогруппах |
| `Orthologs/` | Результаты поиска парных ортологов |
| `Gene_Duplication_Events/Duplications.tsv` | Сведения о событиях дупликации |
| `Species_Tree/` | Эволюционное дерево видов |
| `Resolved_Gene_Trees/` | Разрешённые генетические деревья |

## Распределение белков по техническим категориям

Все **398 017** белков были распределены по техническим категориям на основании результатов OrthoFinder.

| Категория | Количество белков |
|-----------|------------------:|
| `ORTHOLOGS_DETECTED` | 221 956 |
| `ORTHOLOGS_AND_DUPLICATION` | 159 266 |
| `DUPLICATION_DETECTED` | 4 987 |
| `HOMOLOGY_GROUP_ONLY` | 2 182 |
| `UNASSIGNED_CANDIDATE` | 9 626 |
| **Всего** | **398 017** |

Описание категорий:

ORTHOLOGS_DETECTED — белки, для которых найдены записи о парных ортологах.

ORTHOLOGS_AND_DUPLICATION — белки, для которых обнаружены как записи об ортологах, так и события дупликации.

DUPLICATION_DETECTED — белки, связанные с событиями дупликации, но без извлечённой записи о парном ортологе.

HOMOLOGY_GROUP_ONLY — белки, входящие в ортогруппы, но без соответствующей записи о парном ортологе.

UNASSIGNED_CANDIDATE — белки, которые OrthoFinder не смог распределить по ортогруппам.

**Выделение белков *Drosophila melanogaster***

Из общего набора были выделены белки *D. melanogaster*.

Общее количество белковых последовательностей: **13 969** белков.

Их распределение:

| Категория | Количество |
|-----------|-----------:|
| `ORTHOLOGS_DETECTED` | 8 388 |
| `ORTHOLOGS_AND_DUPLICATION` | 5 513 |
| `HOMOLOGY_GROUP_ONLY` | 4 |
| `UNASSIGNED_CANDIDATE` | 64 |
| **Всего** | **13 969** |

Таким образом:

13 901 белок имел сведения об ортологах, включая белки с дополнительной отметкой о дупликации;

4 белка присутствовали в ортогруппах без извлечённых записей о парных ортологах;

64 белка не были распределены по ортогруппам.

Для дальнейшего поиска потенциальных de novo-генов особый интерес представляли последние 64 белка.

**6. Дополнительный поиск гомологии с помощью DIAMOND**

Зачем понадобился дополнительный поиск?

Отсутствие белка в ортогруппах OrthoFinder не обязательно означает отсутствие гомологов у других видов.

Подготовка базы

Из всех протеомов была сформирована единая база, содержащая 398 017 белковых последовательностей.

Для запросов использовались 9 626 белков, не распределённых OrthoFinder по ортогруппам.

Команда DIAMOND

`
diamond makedb \
    --in all_proteomes_db.faa \
    --db all_proteomes_db
`

`
diamond blastp \
    --query unassigned_queries.faa \
    --db all_proteomes_db \
    --ultra-sensitive \
    --evalue 1e-3 \
    --max-target-seqs 100 \
    --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
    --threads 16 \
    --out unassigned_diamond_hits.tsv
`

**Белки *D. melanogaster*, не распределённые по ортогруппам**

Из 64 белков *D. melanogaster*, не распределённых по ортогруппам:

| Результат DIAMOND | Количество |
|-------------------|-----------:|
| `CROSS_SPECIES_DIAMOND_HIT` | 7 |
| `NO_CROSS_SPECIES_DIAMOND_HITS` | 57 |
| **Всего** | **64** |

## Скрипты проекта

Ниже перечислены скрипты, использовавшиеся или подготовленные в ходе работы.

| Скрипт | Назначение | Результаты |
|--------|------------|------------|
| `download_proteomes.py` | Загрузка протеома Dmel из FlyBase и протеомов других видов из NCBI | `proteomes_raw/`, аннотации, `download_report.csv` |
| `classify_gene_origin.py` | Первичная классификация белков по результатам OrthoFinder | `gene_origin_classification.tsv`, `gene_origin_summary.tsv`, `unassigned_candidates.tsv` |
| `prepare_diamond_search.py` | Подготовка FASTA-запросов и объединённой базы для дополнительного поиска | `unassigned_queries.faa`, `all_proteomes_db.faa` |
| `analyze_unassigned_diamond.py` | Обработка результатов DIAMOND и исключение совпадений внутри одного вида | `unassigned_cross_species_hits.tsv`, `unassigned_diamond_summary.tsv`, `unassigned_no_cross_species_hits.txt` |
| `filter_dmel_origin.py` | Выделение белков Dmel и формирование сводок по их классификации и DIAMOND-статусу | `Dmel_all_genes.tsv`, `Dmel_origin_summary.tsv`, списки кандидатов |
| `extract_trio_ids.py` | Извлечение ID белков из ортогрупп трио и рода Drosophila | `trio_all_shared.tsv`, `trio_restricted.tsv`, `drosophila_genus_restricted.tsv`, ID-файлы, `trio_summary.tsv` |
| `extract_dmel_sequences.py` | Извлечение последовательностей Dmel по ID из двух выбранных таблиц | Два FASTA-файла и сводка извлечения |

> **Примечание:** отдельный скрипт подготовки набора `proteomes_1isoform` в имеющейся истории работы однозначно не установлен. Поэтому его имя не указывается. `prepare_dmel_diamond.py` также фигурировал как отдельная подготовка Dmel-специфичных запросов; его не следует включать в основной workflow без проверки, что этот дополнительный поиск действительно запускался.

---

## 11. Итоговая сводка

| Этап | Результат |
|------|-----------|
| Исходный набор всех видов | 398 017 белков |
| Белки *D. melanogaster* | 13 969 |
| Белки Dmel с данными об ортологах | 13 901 |
| Белки Dmel только в ортогруппах | 4 |
| Белки Dmel без ортогрупп | 64 |
| Из 64: с межвидовыми DIAMOND-хитами | 7 |
| Из 64: без межвидовых DIAMOND-хитов | 57 |
| Белки Dmel в ортогруппах всех трёх близких видов | 13 621 |
| Белки Dmel в ортогруппах, ограниченных родом *Drosophila* | 875 |
| Белки Dmel в ортогруппах, ограниченных трио | 19 |

### Основной результат

В результате сравнительно-геномного анализа были выделены **57** белковых последовательностей *D. melanogaster*, которые:

- не были распределены OrthoFinder по ортогруппам;
- не имели межвидовых совпадений в дополнительном поиске DIAMOND.

Эта группа рассматривается как **приоритетный набор потенциальных *de novo*-кандидатов** для дальнейшей проверки.


