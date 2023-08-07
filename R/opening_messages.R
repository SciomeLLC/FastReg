.onAttach <- function(libname, pkgname) {
  msg <- "
  _____                 _     ____                 
 |  ___|   __ _   ___  | |_  |  _ \\    ___    __ _ 
 | |_     / _` | / __| | __| | |_) |  / _ \\  / _` |
 |  _|   | (_| | \\__ \\ | |_  |  _ <  |  __/ | (_| |
 |_|      \\__,_| |___/  \\__| |_| \\_\\  \\___|  \\__, |
                                             |___/

THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"

  # notify the option to update? (5 possible cases)
  # 1. error in getting the current package version -> yes
  # 2. error in getting the github package version -> yes
  # 3. current package version < github package version -> yes
  # 4. current package version > github package version -> no
  # 5. current package version == github package version -> no
  # in short, notify the option to update unless the version numbers match
  # get version of the currently installed package
#   current_pkg_version <- tryCatch(
#     as.character(utils::packageVersion("FastReg")),
#     error = function(e) "unknown")
#   # github url
#   github_url <- paste0(
#     "https://raw.githubusercontent.com/SciomeLLC/FastReg/",
#     "main/DESCRIPTION")
#   # get github description or handle errors
#   github_pkg_desc <- tryCatch(
#     readLines(github_url),
#     warning = function(w) {"github_desc_read_fail"},
#     error = function(e) {"github_desc_read_fail"})
#   # get the version number of github version
#   if (identical(github_pkg_desc, "github_desc_read_fail")) {
#     github_pkg_version <- "unknown"
#   } else {
#     github_pkg_version <- gsub(
#       ".*ersion: ", "", github_pkg_desc[
#         grep("ersion", github_pkg_desc)][1])
#   }
#   # compare versions
#   compare_version_result <- tryCatch(
#     utils::compareVersion(
#       current_pkg_version, github_pkg_version),
#     warning = function(w) {999}, # 999 indicates no need for update
#     error = function(e) {999})
#   # skip update for case 5
#   if (current_pkg_version != "unknown" &
#       github_pkg_version != "unknown" &
#       compare_version_result == 0) {
#     startup_message <- paste0(
#       "Package attached: FastReg v", current_pkg_version,
#       " (same as the most recent version available through GitHub).")
#   } else if (
#     # skip update for case 4
#     current_pkg_version != "unknown" &
#     github_pkg_version != "unknown" &
#     compare_version_result > 0) {
#     startup_message <- paste0(
#       "Package attached: FastReg v", current_pkg_version,
#       " (probably the most recent version available through GitHub).")
#   } else {
#     # simply notify of the OPTION to update the package
#     # this is simply a notification of the option to update,
#     # rather than a recommendation to update
#     startup_message <- paste0(
#       "Package attached: FastReg v", current_pkg_version,
#       "; Most recent version available on GitHub: v", github_pkg_version,
#       "\n\nYou have an option to update the package ",
#       "with the function `update_FastReg()`. ",
#       "If you do so, make sure to restart R.\n")
#   }
  packageStartupMessage(msg)
  if (Sys.info()[['sysname']] == "Darwin" || Sys.info()[['sysname']] == "Windows") {
    mac_warning <- "
**********
This installation of FastReg has detected a Mac which does not support OpenMP.
It will only work in single-threaded mode.
**********
"
    packageStartupMessage(mac_warning)
  }
#   packageStartupMessage(startup_message)
}
