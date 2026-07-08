## Test environments

* Local macOS Sequoia 15.6.1, R 4.6.0
* GitHub Actions macOS latest, R release
* GitHub Actions Ubuntu latest, R release
* GitHub Actions Windows latest, R release

## R CMD check results

0 errors | 0 warnings | 1 note

The note is expected for a first CRAN submission:

* New submission

## Additional checks

* `devtools::test()`: 50 passing tests
* `lintr::lint_package()`: no lints
* `spelling::spell_check_package()`: no spelling errors
