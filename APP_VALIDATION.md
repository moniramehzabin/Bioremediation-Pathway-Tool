# Local app preview validation

This records the original v0.1 interface checks. Updated v0.2 behavior and
validation are documented in APP_UPDATE_v02.md.

The interface was added to repository commit
`f85289e7f33556985fe6444208cd04473a7c5870`.

Completed checks:

- Four automated app checks cover multipart upload parsing, input validation,
  background engine execution, simultaneous-run prevention, failed-domain-job
  warnings, private URL enforcement, cross-origin request rejection, upload
  size limits, sandboxed HTML-report responses and ZIP downloads.
- A full local run of the Micrococcus terreus genome using its saved InterPro
  TSV completed and generated the report through the app backend. It retained
  five supported protocatechuate-branch steps out of six assessed steps.
- The existing repository release audit passed, including protected database
  checks, synthetic annotation controls and existing evidence gates.
- The browser script passed JavaScript syntax validation.

Limitations:

- The Windows launcher and WSL browser-opening behavior have not been executed
  on Windows in the development environment. Test these on the target PC.
- A real browser interaction/visual inspection could not be completed because
  the development browser download was unavailable. UI interactions and layout
  need a target-browser check before release.
- The new full-run check reused existing InterPro results; it did not submit
  another batch of external jobs. The existing web-validation engine is unchanged.
- Initial dependency installation is still required on new computers. This
  preview simplifies running the analysis after setup; it is not a bundled
  installer for arbitrary Windows machines.
