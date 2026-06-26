# Configuration

This ChIP-seq pipeline follows the same configuration pattern as the RNA-seq
pipeline.

Most users should edit only:

```text
config/user_settings.sh
config/metadata.tsv
```

Create them with:

```bash
bash scripts/init_project.sh
```

`config/pipeline_config.sh` is the advanced configuration engine. It loads
`config/user_settings.sh` automatically and defines all derived paths and
defaults used by the scripts.

