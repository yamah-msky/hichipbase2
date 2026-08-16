# hichipbase2

- Author: [`yamah-msky`](https://github.com/yamah-msky)

## Description

An integrated environment for HiChIP analyses, only for myself, mainly powered by Python, pixi.

## HiGlass

Start or recreate the local HiGlass instance with the project-local data and
temporary directories:

```bash
pixi run higlass-start
```

This exposes HiGlass at <http://localhost:8989> and mounts the following
directories into the container:

- `work/results/higlass-data` as `/data`
- `work/results/higlass-tmp` as `/tmp`

Keeping the temporary directory under `work/` is important. `higlass-manage`
uses hard links while staging files, so the input file and
`work/results/higlass-tmp` must be on the same filesystem. Keep inputs under
this repository's `work/` directory, or copy external inputs there before
ingesting them.

For example:

```bash
pixi run higlass-manage ingest \
  "$PWD/work/results/matrix/137_rnap2.mcool"

pixi run higlass-manage list tilesets
```

After a successful ingest and verification, the corresponding staging entry
in `work/results/higlass-tmp` can be removed. The ingested copy is stored under
`work/results/higlass-data`.

Stop and remove the container without deleting the project-local data:

```bash
pixi run higlass-stop
```
