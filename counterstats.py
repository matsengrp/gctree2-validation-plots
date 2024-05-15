from collections import Counter
import random
import statistics


def countermean(c):
    total = sum(k * v for k, v in c.items())
    return total / sum(c.values())

class CounterIndexer:

    def __init__(self, counter: Counter):
        self.counter = counter
        self.items = sorted(counter.items())

    def __getitem__(self, idx):
        counter = 0
        for curr_num, numcount in self.items:
            counter += numcount
            if counter >= idx + 1:
                break
        return curr_num

def countermedian(c):
    indexer = CounterIndexer(c)

    count = sum(c.values())
    half = count // 2
    if count % 2:
        indices = [half]
    else:
        indices = [half - 1, half]
    centervals = [indexer[idx] for idx in indices]
    return sum(centervals) / len(centervals)

def testindexer():
    test_datasets = [
        [random.randint(3, 9) for _ in range(25)],
        [random.randint(3, 9) for _ in range(25)],
        [random.randint(3, 9) for _ in range(24)],
        [random.randint(3, 9) for _ in range(24)],
    ]

    for data in test_datasets:
        dataindexer = CounterIndexer(Counter(data))
        for idx, num in enumerate(sorted(data)):
            if dataindexer[idx] != num:
                print(f"Index {idx} in {data} should be {num} but is {dataindexer[idx]}")

def testmean():
    test_datasets = [
        [random.randint(3, 9) for _ in range(25)],
        [random.randint(3, 9) for _ in range(25)],
        [random.randint(3, 9) for _ in range(24)],
        [random.randint(3, 9) for _ in range(24)],
    ]
    for data in test_datasets:
        test, real = countermean(Counter(data)), statistics.mean(data)
        if test != real:
            print(f"counter mean {test} doesn't match real mean {real} on dataset:\n{sorted(data)}")


def testmedian():
    test_datasets = [
        [random.randint(3, 9) for _ in range(25)],
        [random.randint(3, 9) for _ in range(25)],
        [random.randint(3, 9) for _ in range(24)],
        [random.randint(3, 9) for _ in range(24)],
    ]
    for data in test_datasets:
        test, real = countermedian(Counter(data)), statistics.median(data)
        if test != real:
            print(f"counter median {test} doesn't match real median {real} on dataset:\n{sorted(data)}")

