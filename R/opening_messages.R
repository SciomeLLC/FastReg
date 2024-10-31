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
  packageStartupMessage(msg)
}
