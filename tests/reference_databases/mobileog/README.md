# Tiny test database preparation

This database was built manually from Beatrix 1.6, downloaded from [mobileogdb.flsi.cloud.vt.edu/entries/database_download](https://mobileogdb.flsi.cloud.vt.edu/entries/database_download) (file beatrix-1-6_v1_all.zip). The faa file contained in the zip archive was stripped up to only 2 sequences in mobileOG-db-beatrix-1.6-tiny.faa, and a Diamond database was built using:

```bash
diamond makedb --in mobileOG-db-beatrix-1.6-tiny.faa --db mobileOG_beatrix_1.6_tiny.dmnd
```

Diamond version: 2.1.13
