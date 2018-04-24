package ru.spbau.mit.kazakov.Bioinf;

import org.jetbrains.annotations.NotNull;

import java.util.*;

public class IndexedSet<T> implements Collection<T> {
    private int insertions = 0;
    private Set<T> set = new HashSet<>();
    private Map<T, Integer> map = new HashMap<>();

    @Override
    public int size() {
        return set.size();
    }

    @Override
    public boolean isEmpty() {
        return set.isEmpty();
    }

    @Override
    public boolean contains(Object o) {
        return set.contains(o);
    }

    public int getIndex(@NotNull T o) {
        if (set.contains(o)) {
            return map.get(o);
        } else {
            return -1;
        }
    }

    @NotNull
    @Override
    public Iterator<T> iterator() {
        return set.iterator();
    }

    @NotNull
    @Override
    public Object[] toArray() {
        return set.toArray();
    }

    @NotNull
    @Override
    public <T1> T1[] toArray(@NotNull T1[] t1s) {
        return set.toArray(t1s);
    }

    @Override
    public boolean add(T t) {
        boolean result = set.add(t);
        if(result) {
            map.put(t, insertions);
        }
        insertions++;
        return result;
    }

    @Override
    public boolean remove(Object o) {
        return set.remove(o);
    }

    @Override
    public boolean containsAll(@NotNull Collection<?> collection) {
        return set.containsAll(collection);
    }

    @Override
    public boolean addAll(@NotNull Collection<? extends T> collection) {
        return set.addAll(collection);
    }

    @Override
    public boolean removeAll(@NotNull Collection<?> collection) {
        boolean result = false;
        for(Object o : collection) {
            map.remove(o);
            result |= set.remove(o);
        }

        return result;
    }

    @Override
    public boolean retainAll(@NotNull Collection<?> collection) {
        return false;
    }

    @Override
    public void clear() {
        set.clear();
        map.clear();
        insertions = 0;
    }
}
