package ru.spbau.mit.kazakov.Bioinf;

import org.jetbrains.annotations.NotNull;

import java.util.HashMap;
import java.util.Map;

public enum Activation {
    HCD, CID;

    private static Map<String, Activation> activations = new HashMap<>();
    static {
        activations.put("HCD", Activation.HCD);
        activations.put("CID", Activation.CID);
    }

    @NotNull
    public static Activation getActivationByString(@NotNull String activation) {
        return activations.get(activation);
    }
}
