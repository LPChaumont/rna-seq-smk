import numpy as np
import sys

# modified from https://github.com/Xinglab/rmats-turbo-tutorial/blob/main/scripts/class_exon.py


def get_exon_class(fn):
    if "SE" in fn:
        exon, event_type = exon_SE, "SE"
    elif "RI" in fn:
        exon, event_type = exon_RI, "RI"
    elif "A3SS" in fn:
        exon, event_type = exon_AXSS, "A3SS"
    elif "A5SS" in fn:
        exon, event_type = exon_AXSS, "A5SS"
    elif "MXE" in fn:
        exon, event_type = exon_MXE, "MXE"
    else:
        print(
            "Invalid alternative event type in the input file name. Please modify the file name."
        )
        sys.exit()

    return exon, event_type


class exon_SE(object):
    def __init__(self, line):
        self.line_list = line.replace('"', "").strip().split("\t")
        (
            self.ID,
            self.GeneID,
            self.geneSymbol,
            self.chrom,
            self.strand,
            self.exonStart_0base,
            self.exonEnd,
            self.upstreamES,
            self.upstreamEE,
            self.downstreamES,
            self.downstreamEE,
            self.ID,
            self.IJC_SAMPLE_1,
            self.SJC_SAMPLE_1,
            self.IJC_SAMPLE_2,
            self.SJC_SAMPLE_2,
            self.IncFormLen,
            self.SkipFormLen,
            self.PValue,
            self.FDR,
            self.IncLevel1,
            self.IncLevel2,
            self.IncLevelDifference,
        ) = self.line_list
        self.uniqID = "|".join(
            [
                self.chrom + ":" + self.exonStart_0base + "-" + self.exonEnd,
                self.strand,
                self.upstreamEE,
                self.downstreamES,
            ]
        )
        self.exonStart_0base = int(self.exonStart_0base)
        self.exonEnd = int(self.exonEnd)
        self.upstreamES = int(self.upstreamES)
        self.upstreamEE = int(self.upstreamEE)
        self.downstreamES = int(self.downstreamES)
        self.downstreamEE = int(self.downstreamEE)
        self.IJC_SAMPLE_1 = [int(x) for x in self.IJC_SAMPLE_1.split(",")]
        self.SJC_SAMPLE_1 = [int(x) for x in self.SJC_SAMPLE_1.split(",")]
        self.IJC_SAMPLE_2 = (
            [int(x) for x in self.IJC_SAMPLE_2.split(",")]
            if not self.IJC_SAMPLE_2 == ""
            else []
        )
        self.SJC_SAMPLE_2 = (
            [int(x) for x in self.SJC_SAMPLE_2.split(",")]
            if not self.SJC_SAMPLE_2 == ""
            else []
        )
        self.IncFormLen = int(self.IncFormLen)
        self.SkipFormLen = int(self.SkipFormLen)
        self.PValue = float(self.PValue) if self.PValue != "NA" else "NA"
        self.FDR = float(self.FDR) if self.FDR != "NA" else "NA"
        self.IncLevel1 = (
            [float(x) if x != "NA" else np.nan for x in self.IncLevel1.split(",")]
            if self.IncLevel1 != ""
            else []
        )
        self.IncLevel2 = (
            [float(x) if x != "NA" else np.nan for x in self.IncLevel2.split(",")]
            if self.IncLevel2 != ""
            else []
        )
        self.IncLevelDifference = (
            float(self.IncLevelDifference) if self.IncLevelDifference != "NA" else "NA"
        )
        # calculate average read count
        self.averageCount = float(
            sum(self.IJC_SAMPLE_1)
            + sum(self.SJC_SAMPLE_1)
            + sum(self.IJC_SAMPLE_2)
            + sum(self.SJC_SAMPLE_2)
        ) / float(len(self.IJC_SAMPLE_1) + len(self.IJC_SAMPLE_2))
        self.averageCountSample1 = (
            float(sum(self.IJC_SAMPLE_1) + sum(self.SJC_SAMPLE_1))
            / float(len(self.IJC_SAMPLE_1))
            if not self.IJC_SAMPLE_1 == []
            else np.nan
        )
        self.averageCountSample2 = (
            float(sum(self.IJC_SAMPLE_2) + sum(self.SJC_SAMPLE_2))
            / float(len(self.IJC_SAMPLE_2))
            if not self.IJC_SAMPLE_2 == []
            else np.nan
        )
        self.averageIJC_SAMPLE_1 = (
            float(sum(self.IJC_SAMPLE_1) + sum(self.IJC_SAMPLE_1))
            / float(len(self.IJC_SAMPLE_1))
            if not self.IJC_SAMPLE_1 == []
            else np.nan
        )
        self.averageSJC_SAMPLE_1 = (
            float(sum(self.SJC_SAMPLE_1) + sum(self.SJC_SAMPLE_1))
            / float(len(self.SJC_SAMPLE_1))
            if not self.SJC_SAMPLE_1 == []
            else np.nan
        )
        self.averageIJC_SAMPLE_2 = (
            float(sum(self.IJC_SAMPLE_2) + sum(self.IJC_SAMPLE_2))
            / float(len(self.IJC_SAMPLE_2))
            if not self.IJC_SAMPLE_2 == []
            else np.nan
        )
        self.averageSJC_SAMPLE_2 = (
            float(sum(self.SJC_SAMPLE_2) + sum(self.SJC_SAMPLE_2))
            / float(len(self.SJC_SAMPLE_2))
            if not self.SJC_SAMPLE_2 == []
            else np.nan
        )
        # calculate average PSI
        self.averagePsiSample1 = (
            float(np.nansum(self.IncLevel1)) / len(self.IncLevel1)
            if not self.IncLevel1 == []
            else np.nan
        )
        self.averagePsiSample2 = (
            float(np.nansum(self.IncLevel2)) / len(self.IncLevel2)
            if not self.IncLevel2 == []
            else np.nan
        )
        return

    def __str__(self):
        return (
            "\t".join(
                [self.uniqID] + self.line_list[1:11] + [self.ID] + self.line_list[12:]
            )
            + "\n"
        )


class exon_RI(object):
    def __init__(self, line):
        self.line_list = line.replace('"', "").strip().split("\t")
        (
            self.ID,
            self.GeneID,
            self.geneSymbol,
            self.chrom,
            self.strand,
            self.riExonStart_0base,
            self.riExonEnd,
            self.upstreamES,
            self.upstreamEE,
            self.downstreamES,
            self.downstreamEE,
            self.ID,
            self.IJC_SAMPLE_1,
            self.SJC_SAMPLE_1,
            self.IJC_SAMPLE_2,
            self.SJC_SAMPLE_2,
            self.IncFormLen,
            self.SkipFormLen,
            self.PValue,
            self.FDR,
            self.IncLevel1,
            self.IncLevel2,
            self.IncLevelDifference,
        ) = self.line_list
        self.uniqID = "|".join(
            [
                self.chrom + ":" + self.upstreamEE + "-" + self.downstreamES,
                self.strand,
                self.upstreamES,
                self.downstreamEE,
            ]
        )
        self.riExonStart_0base = int(self.riExonStart_0base)
        self.riExonEnd = int(self.riExonEnd)
        self.upstreamES = int(self.upstreamES)
        self.upstreamEE = int(self.upstreamEE)
        self.downstreamES = int(self.downstreamES)
        self.downstreamEE = int(self.downstreamEE)
        self.IJC_SAMPLE_1 = [int(x) for x in self.IJC_SAMPLE_1.split(",")]
        self.SJC_SAMPLE_1 = [int(x) for x in self.SJC_SAMPLE_1.split(",")]
        self.IJC_SAMPLE_2 = (
            [int(x) for x in self.IJC_SAMPLE_2.split(",")]
            if not self.IJC_SAMPLE_2 == ""
            else []
        )
        self.SJC_SAMPLE_2 = (
            [int(x) for x in self.SJC_SAMPLE_2.split(",")]
            if not self.SJC_SAMPLE_2 == ""
            else []
        )
        self.IncFormLen = int(self.IncFormLen)
        self.SkipFormLen = int(self.SkipFormLen)
        self.PValue = float(self.PValue) if self.PValue != "NA" else "NA"
        self.FDR = float(self.FDR) if self.FDR != "NA" else "NA"
        self.IncLevel1 = (
            [float(x) if x != "NA" else np.nan for x in self.IncLevel1.split(",")]
            if self.IncLevel1 != ""
            else []
        )
        self.IncLevel2 = (
            [float(x) if x != "NA" else np.nan for x in self.IncLevel2.split(",")]
            if self.IncLevel2 != ""
            else []
        )
        self.IncLevelDifference = (
            float(self.IncLevelDifference) if self.IncLevelDifference != "NA" else "NA"
        )
        # calculate average read count
        self.averageCount = float(
            sum(self.IJC_SAMPLE_1)
            + sum(self.SJC_SAMPLE_1)
            + sum(self.IJC_SAMPLE_2)
            + sum(self.SJC_SAMPLE_2)
        ) / float(len(self.IJC_SAMPLE_1) + len(self.IJC_SAMPLE_2))
        self.averageCountSample1 = (
            float(sum(self.IJC_SAMPLE_1) + sum(self.SJC_SAMPLE_1))
            / float(len(self.IJC_SAMPLE_1))
            if not self.IJC_SAMPLE_1 == []
            else np.nan
        )
        self.averageCountSample2 = (
            float(sum(self.IJC_SAMPLE_2) + sum(self.SJC_SAMPLE_2))
            / float(len(self.IJC_SAMPLE_2))
            if not self.IJC_SAMPLE_2 == []
            else np.nan
        )
        self.averageIJC_SAMPLE_1 = (
            float(sum(self.IJC_SAMPLE_1) + sum(self.IJC_SAMPLE_1))
            / float(len(self.IJC_SAMPLE_1))
            if not self.IJC_SAMPLE_1 == []
            else np.nan
        )
        self.averageSJC_SAMPLE_1 = (
            float(sum(self.SJC_SAMPLE_1) + sum(self.SJC_SAMPLE_1))
            / float(len(self.SJC_SAMPLE_1))
            if not self.SJC_SAMPLE_1 == []
            else np.nan
        )
        self.averageIJC_SAMPLE_2 = (
            float(sum(self.IJC_SAMPLE_2) + sum(self.IJC_SAMPLE_2))
            / float(len(self.IJC_SAMPLE_2))
            if not self.IJC_SAMPLE_2 == []
            else np.nan
        )
        self.averageSJC_SAMPLE_2 = (
            float(sum(self.SJC_SAMPLE_2) + sum(self.SJC_SAMPLE_2))
            / float(len(self.SJC_SAMPLE_2))
            if not self.SJC_SAMPLE_2 == []
            else np.nan
        )
        # calculate average PSI
        self.averagePsiSample1 = (
            float(np.nansum(self.IncLevel1)) / len(self.IncLevel1)
            if not self.IncLevel1 == []
            else np.nan
        )
        self.averagePsiSample2 = (
            float(np.nansum(self.IncLevel2)) / len(self.IncLevel2)
            if not self.IncLevel2 == []
            else np.nan
        )
        return

    def __str__(self):
        return (
            "\t".join(
                [self.uniqID] + self.line_list[1:11] + [self.ID] + self.line_list[12:]
            )
            + "\n"
        )


class exon_AXSS(object):
    def __init__(self, line):
        self.line_list = line.replace('"', "").strip().split("\t")
        (
            self.ID,
            self.GeneID,
            self.geneSymbol,
            self.chrom,
            self.strand,
            self.longExonStart_0base,
            self.longExonEnd,
            self.shortES,
            self.shortEE,
            self.flankingES,
            self.flankingEE,
            self.ID,
            self.IJC_SAMPLE_1,
            self.SJC_SAMPLE_1,
            self.IJC_SAMPLE_2,
            self.SJC_SAMPLE_2,
            self.IncFormLen,
            self.SkipFormLen,
            self.PValue,
            self.FDR,
            self.IncLevel1,
            self.IncLevel2,
            self.IncLevelDifference,
        ) = self.line_list

        if int(self.flankingEE) <= int(self.longExonStart_0base):
            self.uniqID = "|".join(
                [
                    self.chrom + ":" + self.flankingES + "-" + self.flankingEE,
                    self.strand,
                    self.longExonStart_0base,
                    self.shortES,
                ]
            )
        elif int(self.flankingES) >= int(self.longExonEnd):
            self.uniqID = "|".join(
                [
                    self.chrom + ":" + self.flankingES + "-" + self.flankingEE,
                    self.strand,
                    self.longExonEnd,
                    self.shortEE,
                ]
            )
        else:
            sys.exit("Error: check A5SS and A3SS file")

        self.longExonStart_0base = int(self.longExonStart_0base)
        self.longExonEnd = int(self.longExonEnd)
        self.shortES = int(self.shortES)
        self.shortEE = int(self.shortEE)
        self.flankingES = int(self.flankingES)
        self.flankingEE = int(self.flankingEE)
        self.IJC_SAMPLE_1 = [int(x) for x in self.IJC_SAMPLE_1.split(",")]
        self.SJC_SAMPLE_1 = [int(x) for x in self.SJC_SAMPLE_1.split(",")]
        self.IJC_SAMPLE_2 = (
            [int(x) for x in self.IJC_SAMPLE_2.split(",")]
            if not self.IJC_SAMPLE_2 == ""
            else []
        )
        self.SJC_SAMPLE_2 = (
            [int(x) for x in self.SJC_SAMPLE_2.split(",")]
            if not self.SJC_SAMPLE_2 == ""
            else []
        )
        self.IncFormLen = int(self.IncFormLen)
        self.SkipFormLen = int(self.SkipFormLen)
        self.PValue = float(self.PValue) if self.PValue != "NA" else "NA"
        self.FDR = float(self.FDR) if self.FDR != "NA" else "NA"
        self.IncLevel1 = (
            [float(x) if x != "NA" else np.nan for x in self.IncLevel1.split(",")]
            if self.IncLevel1 != ""
            else []
        )
        self.IncLevel2 = (
            [float(x) if x != "NA" else np.nan for x in self.IncLevel2.split(",")]
            if self.IncLevel2 != ""
            else []
        )
        self.IncLevelDifference = (
            float(self.IncLevelDifference) if self.IncLevelDifference != "NA" else "NA"
        )
        # calculate average read count
        self.averageCount = float(
            sum(self.IJC_SAMPLE_1)
            + sum(self.SJC_SAMPLE_1)
            + sum(self.IJC_SAMPLE_2)
            + sum(self.SJC_SAMPLE_2)
        ) / float(len(self.IJC_SAMPLE_1) + len(self.IJC_SAMPLE_2))
        self.averageCountSample1 = (
            float(sum(self.IJC_SAMPLE_1) + sum(self.SJC_SAMPLE_1))
            / float(len(self.IJC_SAMPLE_1))
            if not self.IJC_SAMPLE_1 == []
            else np.nan
        )
        self.averageCountSample2 = (
            float(sum(self.IJC_SAMPLE_2) + sum(self.SJC_SAMPLE_2))
            / float(len(self.IJC_SAMPLE_2))
            if not self.IJC_SAMPLE_2 == []
            else np.nan
        )
        self.averageIJC_SAMPLE_1 = (
            float(sum(self.IJC_SAMPLE_1) + sum(self.IJC_SAMPLE_1))
            / float(len(self.IJC_SAMPLE_1))
            if not self.IJC_SAMPLE_1 == []
            else np.nan
        )
        self.averageSJC_SAMPLE_1 = (
            float(sum(self.SJC_SAMPLE_1) + sum(self.SJC_SAMPLE_1))
            / float(len(self.SJC_SAMPLE_1))
            if not self.SJC_SAMPLE_1 == []
            else np.nan
        )
        self.averageIJC_SAMPLE_2 = (
            float(sum(self.IJC_SAMPLE_2) + sum(self.IJC_SAMPLE_2))
            / float(len(self.IJC_SAMPLE_2))
            if not self.IJC_SAMPLE_2 == []
            else np.nan
        )
        self.averageSJC_SAMPLE_2 = (
            float(sum(self.SJC_SAMPLE_2) + sum(self.SJC_SAMPLE_2))
            / float(len(self.SJC_SAMPLE_2))
            if not self.SJC_SAMPLE_2 == []
            else np.nan
        )
        # calculate average PSI
        self.averagePsiSample1 = (
            float(np.nansum(self.IncLevel1)) / len(self.IncLevel1)
            if not self.IncLevel1 == []
            else np.nan
        )
        self.averagePsiSample2 = (
            float(np.nansum(self.IncLevel2)) / len(self.IncLevel2)
            if not self.IncLevel2 == []
            else np.nan
        )
        return

    def __str__(self):
        return (
            "\t".join(
                [self.uniqID] + self.line_list[1:11] + [self.ID] + self.line_list[12:]
            )
            + "\n"
        )


class exon_MXE(object):
    def __init__(self, line):
        self.line_list = line.replace('"', "").strip().split("\t")
        (
            self.ID,
            self.GeneID,
            self.geneSymbol,
            self.chrom,
            self.strand,
            self.stExonStart_0base,
            self.stExonEnd,
            self.ndExonStart_0base,
            self.ndExonEnd,
            self.upstreamES,
            self.upstreamEE,
            self.downstreamES,
            self.downstreamEE,
            self.ID,
            self.IJC_SAMPLE_1,
            self.SJC_SAMPLE_1,
            self.IJC_SAMPLE_2,
            self.SJC_SAMPLE_2,
            self.IncFormLen,
            self.SkipFormLen,
            self.PValue,
            self.FDR,
            self.IncLevel1,
            self.IncLevel2,
            self.IncLevelDifference,
        ) = self.line_list

        self.uniqID = "|".join(
            [
                self.chrom
                + ":"
                + self.stExonStart_0base
                + "-"
                + self.stExonEnd
                + ":"
                + self.ndExonStart_0base
                + "-"
                + self.ndExonEnd,
                self.strand,
                self.upstreamEE,
                self.downstreamES,
            ]
        )
        self.stExonStart_0base = int(self.stExonStart_0base)
        self.stExonEnd = int(self.stExonEnd)
        self.ndExonStart_0base = int(self.ndExonStart_0base)
        self.ndExonEnd = int(self.ndExonEnd)
        self.upstreamES = int(self.upstreamES)
        self.upstreamEE = int(self.upstreamEE)
        self.downstreamES = int(self.downstreamES)
        self.downstreamEE = int(self.downstreamEE)
        self.IJC_SAMPLE_1 = [int(x) for x in self.IJC_SAMPLE_1.split(",")]
        self.SJC_SAMPLE_1 = [int(x) for x in self.SJC_SAMPLE_1.split(",")]
        self.IJC_SAMPLE_2 = (
            [int(x) for x in self.IJC_SAMPLE_2.split(",")]
            if not self.IJC_SAMPLE_2 == ""
            else []
        )
        self.SJC_SAMPLE_2 = (
            [int(x) for x in self.SJC_SAMPLE_2.split(",")]
            if not self.SJC_SAMPLE_2 == ""
            else []
        )
        self.IncFormLen = int(self.IncFormLen)
        self.SkipFormLen = int(self.SkipFormLen)
        self.PValue = float(self.PValue) if self.PValue != "NA" else "NA"
        self.FDR = float(self.FDR) if self.FDR != "NA" else "NA"
        self.IncLevel1 = (
            [float(x) if x != "NA" else np.nan for x in self.IncLevel1.split(",")]
            if self.IncLevel1 != ""
            else []
        )
        self.IncLevel2 = (
            [float(x) if x != "NA" else np.nan for x in self.IncLevel2.split(",")]
            if self.IncLevel2 != ""
            else []
        )
        self.IncLevelDifference = (
            float(self.IncLevelDifference) if self.IncLevelDifference != "NA" else "NA"
        )
        # calculate average read count
        self.averageCount = float(
            sum(self.IJC_SAMPLE_1)
            + sum(self.SJC_SAMPLE_1)
            + sum(self.IJC_SAMPLE_2)
            + sum(self.SJC_SAMPLE_2)
        ) / float(len(self.IJC_SAMPLE_1) + len(self.IJC_SAMPLE_2))
        self.averageCountSample1 = (
            float(sum(self.IJC_SAMPLE_1) + sum(self.SJC_SAMPLE_1))
            / float(len(self.IJC_SAMPLE_1))
            if not self.IJC_SAMPLE_1 == []
            else np.nan
        )
        self.averageCountSample2 = (
            float(sum(self.IJC_SAMPLE_2) + sum(self.SJC_SAMPLE_2))
            / float(len(self.IJC_SAMPLE_2))
            if not self.IJC_SAMPLE_2 == []
            else np.nan
        )
        self.averageIJC_SAMPLE_1 = (
            float(sum(self.IJC_SAMPLE_1) + sum(self.IJC_SAMPLE_1))
            / float(len(self.IJC_SAMPLE_1))
            if not self.IJC_SAMPLE_1 == []
            else np.nan
        )
        self.averageSJC_SAMPLE_1 = (
            float(sum(self.SJC_SAMPLE_1) + sum(self.SJC_SAMPLE_1))
            / float(len(self.SJC_SAMPLE_1))
            if not self.SJC_SAMPLE_1 == []
            else np.nan
        )
        self.averageIJC_SAMPLE_2 = (
            float(sum(self.IJC_SAMPLE_2) + sum(self.IJC_SAMPLE_2))
            / float(len(self.IJC_SAMPLE_2))
            if not self.IJC_SAMPLE_2 == []
            else np.nan
        )
        self.averageSJC_SAMPLE_2 = (
            float(sum(self.SJC_SAMPLE_2) + sum(self.SJC_SAMPLE_2))
            / float(len(self.SJC_SAMPLE_2))
            if not self.SJC_SAMPLE_2 == []
            else np.nan
        )
        # calculate average PSI
        self.averagePsiSample1 = (
            float(np.nansum(self.IncLevel1)) / len(self.IncLevel1)
            if not self.IncLevel1 == []
            else np.nan
        )
        self.averagePsiSample2 = (
            float(np.nansum(self.IncLevel2)) / len(self.IncLevel2)
            if not self.IncLevel2 == []
            else np.nan
        )
        return

    def __str__(self):
        return (
            "\t".join(
                [self.uniqID] + self.line_list[1:13] + [self.ID] + self.line_list[14:]
            )
            + "\n"
        )
