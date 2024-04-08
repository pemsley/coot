
dir=$1

for file in $dir/* ;
do
    if test -d $file ; then
        :
    else
        sed -ibackup -e 's/it under the terms of the GNU General Public License as published by/it under the terms of the GNU Lesser General Public License as published/' $file

        sed -ibackup-2 -e 's/\* the Free Software Foundation; either version 3 of the License, or/* by the Free Software Foundation; either version 3 of the License, or/' $file

        sed -ibackup-3 -e 's/\* General Public License for more details/* Lesser General Public License for more details/' $file

        sed -ibackup-4 -e 's/You should have received a copy of the GNU General Public License$/& and/' $file
        sed -ibackup-5 -e 's/along with this program; if not, write to the Free Software/the GNU Lesser General Public License along with this program; if not,/' $file
        sed -ibackup-6 -e 's/Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA/write to the Free Software Foundation, Inc., 51 Franklin Street,/' $file
        sed -ibackup-7 -e 's/\* 02110-1301, USA/* Fifth Floor, Boston, MA, 02110-1301, USA./' $file
    fi
done


 # * You should have received a copy of the GNU General Public License
 # * along with this program; if not, write to the Free Software
 # * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 # * 02110-1301, USA

 # * You should have received a copies of the GNU General Public License and
 # * the GNU Lesser General Public License along with this program; if not,
 # * write to the Free Software Foundation, Inc., 51 Franklin Street,
 # * Fifth Floor, Boston, MA, 02110-1301, USA.
 # * See http://www.gnu.org/licenses/
 # *
