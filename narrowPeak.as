table narrowPeak
"BED6+4 Peaks of signal enrichment"
(
    string chrom;       "Reference sequence chromosome or scaffold"
    uint   chromStart;  "Start position in chromosome"
    uint   chromEnd;    "End position in chromosome"
    string name;        "Name given to a region"
    uint   score;       "Score 0-1000"
    char[1] strand;     "+ or - or ."
    float  signalValue; "Measurement of average enrichment"
    float  pValue;      "Significance -log10"
    float  qValue;      "Significance -log10"
    int    peak;        "Summit position"
)
