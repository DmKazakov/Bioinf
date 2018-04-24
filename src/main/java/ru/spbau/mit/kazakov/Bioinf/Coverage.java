package ru.spbau.mit.kazakov.Bioinf;

import org.jetbrains.annotations.NotNull;

public class Coverage {
    private String hcdHeavyChain;
    private String cidHeavyChain;
    private String hcdLightChain;
    private String cidLightChain;

    public Coverage(@NotNull String hcdHeavyChain, @NotNull String cidHeavyChain, @NotNull String hcdLightChain, @NotNull String cidLightChain) {
        this.hcdHeavyChain = hcdHeavyChain;
        this.cidHeavyChain = cidHeavyChain;
        this.hcdLightChain = hcdLightChain;
        this.cidLightChain = cidLightChain;
    }

    @NotNull
    public String getHcdHeavyChain() {
        return hcdHeavyChain;
    }

    @NotNull
    public String getCidHeavyChain() {
        return cidHeavyChain;
    }

    @NotNull
    public String getHcdLightChain() {
        return hcdLightChain;
    }

    @NotNull
    public String getCidLightChain() {
        return cidLightChain;
    }
}
