package ru.spbau.mit.kazakov.Bioinf;

import org.jetbrains.annotations.NotNull;

public class Sequence {
    private String sequence;
    private double rnkScr;

    public Sequence(@NotNull String sequence, double rnkScr) {
        this.sequence = sequence;
        this.rnkScr = rnkScr;
    }

    @NotNull
    public String getSequence() {
        return sequence;
    }

    public double getRnkScr() {
        return rnkScr;
    }
}
