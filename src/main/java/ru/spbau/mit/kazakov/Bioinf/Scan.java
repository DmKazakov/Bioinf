package ru.spbau.mit.kazakov.Bioinf;

import org.jetbrains.annotations.NotNull;

import java.util.List;

public class Scan {
    private int number;
    private Activation activation;
    private List<Sequence> sequences;
    private Scan pair = null;
    private double precursorMass;

    public Scan(int number, @NotNull List<Sequence> sequences) {
        this.number = number;
        this.sequences = sequences;
    }

    @NotNull
    public List<Sequence> getSequences() {
        return sequences;
    }

    @NotNull
    public Activation getActivation() {
        return activation;
    }

    public void setActivation(@NotNull Activation activation) {
        this.activation = activation;
    }

    public int getNumber() {
        return number;
    }

    public boolean hasPair() {
        return pair == null;
    }

    @NotNull
    public Scan getPair() {
        return pair;
    }

    public double getPrecursorMass() {
        return precursorMass;
    }

    public void setPrecursorMass(double precursorMass) {
        this.precursorMass = precursorMass;
    }
}
