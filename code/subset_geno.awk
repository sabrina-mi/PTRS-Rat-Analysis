BEGIN { OFS="\t" }
NR==FNR {
    samples[++numSamples] = $1
    next
}
FNR==1 {
    for (i=1; i<=NF; i++) {
        f[$i] = i
    }
}
{
    printf "%s%s%s%s", $1, OFS, $2, OFS
    for (sampleNr =1; sampleNr <= nSamples; sampleNr++) {
        sample = samples[sampleNr]
        printf "%s%s", $(f[sample]), (sampleNr<numSamples ? OFS : ORS)
    }
}