# Contributing

## Regenerating documentation

The help pages (`man/*.Rd`) and `NAMESPACE` are generated from the roxygen
comments in `R/`. Do not edit them by hand. Use roxygen2 8.1.0 or newer (it
needs R >= 4.1); the version is recorded in `DESCRIPTION` as
`Config/roxygen2/version`. From the repository root:

```r
roxygen2::roxygenise()      # or devtools::document()
```

Commit the changes to `man/` and `NAMESPACE` together with the code. Older
roxygen2 versions (7.x) will rewrite the version field and warn about the
`@exportS3Method NULL` tags, so check `packageVersion("roxygen2")` first.