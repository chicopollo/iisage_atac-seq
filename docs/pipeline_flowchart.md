# IISAGE ATAC-seq Pipeline v2.0 - Flowchart

## Complete Pipeline Architecture

```mermaid
flowchart TD
    Start([Pipeline Start]) --> CheckMode{Deployment<br/>Mode?}

    %% Deployment Branch
    CheckMode -->|DEPLOY_ALL=true| FindSpecies[Find All Species<br/>with config.yaml]
    FindSpecies --> DeployJobs[Submit SLURM Jobs<br/>for Each Species]
    DeployJobs --> EndDeploy([Deployment Complete])

    %% Single Species Branch
    CheckMode -->|Single Species| LoadConfig[Load Species Config<br/>config.yaml]
    LoadConfig --> ParseYAML[Safe YAML Parsing<br/>No eval]
    ParseYAML --> ValidateConfig{Valid Config?}
    ValidateConfig -->|No| Error1([Exit: Missing Config])
    ValidateConfig -->|Yes| CheckDisk[Check Disk Space<br/>≥100GB free]

    CheckDisk --> LoadModules[Load & Validate<br/>Bioinformatics Modules]
    LoadModules --> ModulesOK{Modules<br/>Loaded?}
    ModulesOK -->|No| Error2([Exit: Module Failed])
    ModulesOK -->|Yes| LoadR[Load R Environment<br/>Optional]

    LoadR --> ValidatePaths[Validate Paths<br/>Samples, Reference]
    ValidatePaths --> CreateDirs[Create Output<br/>Directories]
    CreateDirs --> RecordVersions[Record Software<br/>Versions]
    RecordVersions --> CheckArray{Array Job?}

    %% Array vs Sequential Processing
    CheckArray -->|Yes| ArrayMode[Array Task Mode<br/>SLURM_ARRAY_TASK_ID set]
    CheckArray -->|No| SeqMode[Sequential Mode<br/>Main Job]

    %% Stage 1: Trimming & Alignment
    ArrayMode --> Stage1A[STAGE 1: Trimming & Alignment]
    SeqMode --> Stage1S[STAGE 1: Trimming & Alignment]

    Stage1A --> BuildIndex1[Build/Validate<br/>Bowtie2 Index]
    Stage1S --> BuildIndex2[Build/Validate<br/>Bowtie2 Index]

    BuildIndex1 --> GetSamples1[Get Sample List]
    BuildIndex2 --> GetSamples2[Get Sample List]

    GetSamples1 --> ProcessArraySample[Process Single Sample<br/>Task ID → Sample]
    GetSamples2 --> ProcessAllSamples[Process All Samples<br/>Sequential]

    ProcessArraySample --> MergeFASTQ1[Merge FASTQ Files<br/>zcat *.fq.gz]
    ProcessAllSamples --> MergeFASTQ2[Merge FASTQ Files<br/>zcat *.fq.gz]

    MergeFASTQ1 --> Trim1[Trimmomatic<br/>Adapter Removal]
    MergeFASTQ2 --> Trim2[Trimmomatic<br/>Adapter Removal]

    Trim1 --> Align1[Bowtie2 Alignment<br/>--very-sensitive]
    Trim2 --> Align2[Bowtie2 Alignment<br/>--very-sensitive]

    Align1 --> SAMtoBAM1[SAMtools<br/>Convert & Sort BAM]
    Align2 --> SAMtoBAM2[SAMtools<br/>Convert & Sort BAM]

    SAMtoBAM1 --> IndexBAM1[SAMtools Index<br/>Create .bai]
    SAMtoBAM2 --> IndexBAM2[SAMtools Index<br/>Create .bai]

    IndexBAM1 --> ValidateBAM1[Validate BAM<br/>samtools quickcheck]
    IndexBAM2 --> ValidateBAM2[Validate BAM<br/>samtools quickcheck]

    ValidateBAM1 --> Checkpoint1[Create Checkpoint<br/>sample_aligned.done]
    ValidateBAM2 --> Checkpoint2[Create Checkpoint<br/>sample_aligned.done]

    Checkpoint1 --> ArrayExit{Array Task?}
    Checkpoint2 --> Stage2Check{Run Peak<br/>Calling?}

    ArrayExit -->|Yes| ExitArray([Array Task Complete])
    ArrayExit -->|No| Stage2Check

    %% Stage 2: Peak Calling
    Stage2Check -->|No| Stage3Check{Run<br/>DeepTools?}
    Stage2Check -->|Yes| Stage2[STAGE 2: Peak Calling]

    Stage2 --> IterateSamples[For Each Sample<br/>BAM File]
    IterateSamples --> SizeSelect[Size Selection<br/>Filter ≤100bp fragments]
    SizeSelect --> FilterChrom[Filter Chromosomes<br/>Remove random, chrUn, chrM]
    FilterChrom --> FilterMAPQ[Filter Quality<br/>MAPQ ≥ 5]
    FilterMAPQ --> SortBAM[Sort & Index<br/>Size-selected BAM]

    SortBAM --> MACS2[MACS2 Peak Calling<br/>--nomodel --call-summits]
    MACS2 --> PeakFiles[Generate Outputs<br/>narrowPeak, summits, xls]
    PeakFiles --> Checkpoint3[Create Checkpoint<br/>sample_peaks.done]
    Checkpoint3 --> MoreSamples{More<br/>Samples?}
    MoreSamples -->|Yes| IterateSamples
    MoreSamples -->|No| Stage3Check

    %% Stage 3: DeepTools
    Stage3Check -->|No| Stage4Check{Run R<br/>Analysis?}
    Stage3Check -->|Yes| Stage3[STAGE 3: DeepTools]

    Stage3 --> CollectBAMs[Collect All<br/>Size-selected BAMs]
    CollectBAMs --> CheckTSS{TSS BED<br/>Available?}

    CheckTSS -->|Yes| ComputeMatrix[computeMatrix<br/>TSS ±3kb regions]
    ComputeMatrix --> PlotHeatmap[plotHeatmap<br/>Visualize TSS]
    PlotHeatmap --> MultiBam[multiBamSummary<br/>Genome-wide bins]

    CheckTSS -->|No| MultiBam

    MultiBam --> PlotCorr[plotCorrelation<br/>Pearson correlation]
    PlotCorr --> CreateBigWig[bamCoverage<br/>Generate BigWig files]
    CreateBigWig --> Checkpoint4[Create Checkpoint<br/>deeptools_complete.done]
    Checkpoint4 --> Stage4Check

    %% Stage 4: R Analysis
    Stage4Check -->|No| Summary[Generate Pipeline<br/>Summary]
    Stage4Check -->|Yes| Stage4[STAGE 4: R Analysis]

    Stage4 --> CheckPeaks{Peaks &<br/>BAMs exist?}
    CheckPeaks -->|No| SkipR[Skip R Analysis<br/>Log Warning]
    CheckPeaks -->|Yes| ValidateMeta[Validate<br/>sample_metadata.csv]

    ValidateMeta --> MetaOK{Metadata<br/>Valid?}
    MetaOK -->|No| MetaRequired{R Analysis<br/>Required?}
    MetaRequired -->|Yes| Error3([Exit: No Metadata])
    MetaRequired -->|No| SkipR

    MetaOK -->|Yes| RunR[Execute R Script<br/>generate_atac_figs_v2.R]

    RunR --> RVenn[Generate Venn Diagrams<br/>All & per-condition]
    RVenn --> RCorr[Correlation Heatmaps<br/>DiffBind analysis]
    RCorr --> RAnnot{Skip<br/>Annotation?}

    RAnnot -->|No| LoadGTF[Load GTF/GFF<br/>Annotation]
    LoadGTF --> AnnotatePeaks[ChIPseeker<br/>Peak Annotation]
    AnnotatePeaks --> GenomicDist[Genomic Distribution<br/>Promoter/Exon/Intron]
    GenomicDist --> TSSdist[Distance to TSS<br/>Analysis]
    TSSdist --> RMotif{Skip<br/>Motifs?}

    RAnnot -->|Yes| RMotif

    RMotif -->|No| HOMER[HOMER Motif Analysis<br/>findMotifsGenome.pl]
    HOMER --> MotifHTML[Generate HTML Reports<br/>Known & de novo motifs]
    MotifHTML --> RFeatures[Peak Feature Tables<br/>Width, signal, counts]

    RMotif -->|Yes| RFeatures

    RFeatures --> Consensus[Consensus Peak Sets<br/>All & per-condition]
    Consensus --> ROK{R Analysis<br/>Success?}

    ROK -->|Yes| Checkpoint5[Create Checkpoint<br/>r_analysis_complete.done]
    Checkpoint5 --> EmailR[Send Success Email<br/>List all outputs]
    EmailR --> Summary

    ROK -->|No| RReq{R Analysis<br/>Required?}
    RReq -->|Yes| Error4([Exit: R Failed])
    RReq -->|No| EmailRFail[Send Failure Email<br/>Log error location]

    SkipR --> Summary
    EmailRFail --> Summary

    %% Summary & Completion
    Summary --> CountReads[Count Aligned Reads<br/>Per sample statistics]
    CountReads --> CountPeaks[Count Called Peaks<br/>Per sample statistics]
    CountPeaks --> CountFigures{R Analysis<br/>Completed?}

    CountFigures -->|Yes| CountROutputs[Count R Outputs<br/>Figures, tables, motifs]
    CountFigures -->|No| WriteSummary[Write Summary File<br/>pipeline_summary.txt]
    CountROutputs --> WriteSummary

    WriteSummary --> LogDirs[Log Output<br/>Directories]
    LogDirs --> FinalEmail[Send Completion Email<br/>Include summary]
    FinalEmail --> End([Pipeline Complete])

    %% Styling
    classDef stageClass fill:#e1f5ff,stroke:#01579b,stroke-width:3px
    classDef decisionClass fill:#fff9c4,stroke:#f57f17,stroke-width:2px
    classDef errorClass fill:#ffcdd2,stroke:#c62828,stroke-width:2px
    classDef processClass fill:#c8e6c9,stroke:#2e7d32,stroke-width:2px
    classDef checkpointClass fill:#e1bee7,stroke:#6a1b9a,stroke-width:2px

    class Stage1A,Stage1S,Stage2,Stage3,Stage4 stageClass
    class CheckMode,ValidateConfig,ModulesOK,CheckArray,Stage2Check,Stage3Check,Stage4Check,CheckPeaks,MetaOK,MetaRequired,CheckTSS,RAnnot,RMotif,ROK,RReq,CountFigures,ArrayExit,MoreSamples decisionClass
    class Error1,Error2,Error3,Error4 errorClass
    class Checkpoint1,Checkpoint2,Checkpoint3,Checkpoint4,Checkpoint5 checkpointClass
```

## Simplified High-Level Flow

```mermaid
graph LR
    A[Raw FASTQ<br/>Files] -->|Trimmomatic| B[Trimmed<br/>Reads]
    B -->|Bowtie2| C[Aligned<br/>BAM]
    C -->|Filter & Sort| D[Processed<br/>BAM]
    D -->|MACS2| E[Peak<br/>Calls]
    E -->|deepTools| F[Visualizations<br/>BigWig]
    E -->|R Analysis| G[Figures &<br/>Tables]
    D -->|R Analysis| G

    style A fill:#ffebee
    style B fill:#e3f2fd
    style C fill:#e8f5e9
    style D fill:#f3e5f5
    style E fill:#fff3e0
    style F fill:#e0f2f1
    style G fill:#fce4ec
```

## Data Flow by Stage

```mermaid
graph TD
    subgraph "Input Data"
        FASTQ[FASTQ Files<br/>*_1.fq.gz, *_2.fq.gz]
        CONFIG[config.yaml<br/>Species Configuration]
        META[sample_metadata.csv<br/>Sample Information]
        REF[Reference Genome<br/>FASTA + GTF]
    end

    subgraph "Stage 1: Trimming & Alignment"
        TRIM[Trimmed Reads<br/>paired.fastq]
        SAM[Aligned Reads<br/>aligned.sam]
        BAM[Sorted BAM<br/>sorted.bam + .bai]
    end

    subgraph "Stage 2: Peak Calling"
        FILT[Size-selected BAM<br/>≤100bp fragments]
        PEAKS[Peak Files<br/>narrowPeak, xls, summits]
    end

    subgraph "Stage 3: DeepTools"
        BIGWIG[BigWig Tracks<br/>.bw files]
        MATRIX[Coverage Matrix<br/>TSS regions]
        CORR[Correlation Plots<br/>PNG]
    end

    subgraph "Stage 4: R Analysis"
        VENN[Venn Diagrams<br/>PNG]
        ANNOT[Peak Annotations<br/>CSV]
        MOTIFS[Motif Results<br/>HTML]
        TABLES[Feature Tables<br/>CSV]
        CONSENSUS[Consensus Peaks<br/>BED]
    end

    subgraph "Final Outputs"
        SUMMARY[Pipeline Summary<br/>TXT]
        VERSIONS[Software Versions<br/>TXT]
        FIGS[Figures Directory<br/>All visualizations]
    end

    FASTQ --> TRIM
    CONFIG --> TRIM
    TRIM --> SAM
    SAM --> BAM
    BAM --> FILT
    FILT --> PEAKS
    FILT --> BIGWIG
    FILT --> MATRIX
    MATRIX --> CORR
    PEAKS --> VENN
    PEAKS --> ANNOT
    META --> VENN
    REF --> ANNOT
    PEAKS --> MOTIFS
    REF --> MOTIFS
    PEAKS --> TABLES
    PEAKS --> CONSENSUS
    VENN --> FIGS
    ANNOT --> FIGS
    MOTIFS --> FIGS
    TABLES --> FIGS
    CONSENSUS --> FIGS
    BAM --> SUMMARY
    PEAKS --> SUMMARY
    FIGS --> SUMMARY
    CONFIG --> VERSIONS
```

## Parallel Processing (Array Jobs)

```mermaid
graph TD
    Main[Main Job<br/>SLURM Submission] --> Array[Array Job Launched<br/>--array=1-N]

    Array --> Task1[Task 1<br/>Process Sample 1]
    Array --> Task2[Task 2<br/>Process Sample 2]
    Array --> TaskN[Task N<br/>Process Sample N]

    Task1 --> T1Done[Checkpoint<br/>sample1_aligned.done]
    Task2 --> T2Done[Checkpoint<br/>sample2_aligned.done]
    TaskN --> TNDone[Checkpoint<br/>sampleN_aligned.done]

    T1Done --> Wait{All Tasks<br/>Complete?}
    T2Done --> Wait
    TNDone --> Wait

    Wait -->|Yes| MainResumes[Main Job Resumes<br/>No SLURM_ARRAY_TASK_ID]
    MainResumes --> PeakCalling[Stage 2: Peak Calling<br/>All samples]
    PeakCalling --> DeepTools[Stage 3: DeepTools<br/>All samples]
    DeepTools --> RAnalysis[Stage 4: R Analysis<br/>All samples together]
    RAnalysis --> Complete[Pipeline Complete]

    style Task1 fill:#bbdefb
    style Task2 fill:#c5cae9
    style TaskN fill:#d1c4e9
    style MainResumes fill:#c8e6c9
    style RAnalysis fill:#ffe0b2
```

## Checkpoint System

```mermaid
stateDiagram-v2
    [*] --> CheckResume: Pipeline Start

    CheckResume --> Stage1: No checkpoint/<br/>Resume disabled
    CheckResume --> SkipStage1: Checkpoint exists

    Stage1 --> CreateCP1: Stage complete
    CreateCP1 --> Stage2

    SkipStage1 --> Stage2

    Stage2 --> CheckCP2: Check sample_peaks.done
    CheckCP2 --> RunStage2: Not found
    CheckCP2 --> SkipStage2: Found

    RunStage2 --> CreateCP2: Stage complete
    CreateCP2 --> Stage3

    SkipStage2 --> Stage3

    Stage3 --> CheckCP3: Check deeptools_complete.done
    CheckCP3 --> RunStage3: Not found
    CheckCP3 --> SkipStage3: Found

    RunStage3 --> CreateCP3: Stage complete
    CreateCP3 --> Stage4

    SkipStage3 --> Stage4

    Stage4 --> CheckCP4: Check r_analysis_complete.done
    CheckCP4 --> RunStage4: Not found
    CheckCP4 --> SkipStage4: Found

    RunStage4 --> CreateCP4: Stage complete
    CreateCP4 --> [*]

    SkipStage4 --> [*]
```

## Error Handling Flow

```mermaid
graph TD
    Process[Pipeline Process] --> Error{Error<br/>Occurs?}

    Error -->|No| Continue[Continue Pipeline]

    Error -->|Yes| ErrorType{Error<br/>Type?}

    ErrorType -->|Critical| CheckRequired{Stage<br/>Required?}
    ErrorType -->|Warning| LogWarn[Log Warning]

    CheckRequired -->|Yes| SendFail[Send Failure Email]
    CheckRequired -->|No| LogError[Log Error]

    SendFail --> Exit[Exit Pipeline<br/>Status: Failed]

    LogWarn --> Continue
    LogError --> Continue

    Continue --> NextStage[Next Stage]

    style Error fill:#fff9c4
    style ErrorType fill:#ffe0b2
    style SendFail fill:#ffcdd2
    style Exit fill:#ef9a9a
    style LogWarn fill:#fff59d
```

---

## Legend

| Symbol | Meaning |
|--------|---------|
| 🟦 **Stage** | Main pipeline stage (1-4) |
| 🟨 **Decision** | Conditional branching point |
| 🟥 **Error** | Error exit point |
| 🟩 **Process** | Data processing step |
| 🟪 **Checkpoint** | Resume capability marker |

## Key Features Highlighted

1. **Array Processing:** Samples processed in parallel for speed
2. **Checkpoint System:** Resume from failure without reprocessing
3. **Modular Stages:** Each stage can be enabled/disabled independently
4. **Error Resilience:** Non-critical failures don't stop pipeline
5. **Validation:** Multiple validation steps ensure data integrity

