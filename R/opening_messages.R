.onAttach <- function(libname, pkgname) {
  msg <- "
,------.               ,--.   
|  .---',--,--. ,---.,-'  '-. 
|  `--,' ,-.  |(  .-''-.  .-' 
|  |`  \\ '-'  |.-'  `) |  |   
`--'    `--`--'`----'  `--'   
,--.   ,--.,--.     ,---.     
 \\  `.'  / |  |    /  O  \\    
  \\     /  |  |   |  .-.  |   
   \\   /   |  '--.|  | |  |   
    `-'    `-----'`--' `--'   

THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"
  packageStartupMessage(msg)
  if (Sys.info()[["sysname"]] == "Darwin") {
    mac_warning <- "
**********
This installation of FastReg has detected a Mac which does not support OpenMP.
It will only work in single-threaded mode.
**********
"
    packageStartupMessage(mac_warning)
  }
}
