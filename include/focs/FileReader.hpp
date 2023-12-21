
#ifndef FOCS_FILEREADER
#define FOCS_FILEREADER

#include <array>
#include <cstdio>
#include <string>
#include <unordered_map>

#include <boost/filesystem.hpp>

#include "DataProvider.hpp"
#include "DataRecord.hpp"

namespace focs {

extern const std::array<int, 4> NETCDF1_MAGIC_NUMBER;
extern const std::array<int, 4> NETCDF2_MAGIC_NUMBER;
extern const std::array<int, 4> HDF4_MAGIC_NUMBER;
extern const std::array<int, 8> HDF5_MAGIC_NUMBER;

// Everything that writes to hints.attributes_ is performing a copy (of what they're adding).  Switch to smart pointers?
class FileReaderHints {
    public:
        FileReaderHints() = default;
        explicit FileReaderHints(const std::string& path);

        bool magic_number(const std::string& bytes) const noexcept;
        template<std::size_t N> bool magic_number(const std::array<int, N>& bytes) const noexcept;

        const std::array<int, 8>& first_bytes() const noexcept {return first_bytes_;}
        void first_bytes(const std::array<int, 8>& bytes){first_bytes_ = bytes;}

        uintmax_t size() const noexcept {return size_;}
        void size(uintmax_t size){size_ = size;}

        const std::string& filename() const noexcept {return filename_;}
        void filename(const std::string& filename){filename_ = filename;}

        const std::string& format() const noexcept {return format_;}
        void format(const std::string& format){format_ = format;}

        auto& attributes(){return attributes_;}

        bool is_hdf4() const noexcept;
        bool is_hdf5() const noexcept;
        bool is_netcdf() const noexcept;
    private:
        std::string filename_{};
        std::array<int, 8> first_bytes_{{0, 0, 0, 0, 0, 0, 0, 0}};
        uintmax_t size_{0};
        std::string format_{"unknown"};
        std::unordered_map<std::string, std::string> attributes_{};
};
class FileReader : public DataProvider {
    public:
        FileReader(const std::string& name, const std::string& description) : DataProvider(name, description) {}
        FileReader(const std::string& name) : FileReader(name, {}) {}

        enum validity {
            invalid, possibly_valid, valid
        };
        virtual ~FileReader() override {}

        // virtual const std::vector<Product>& provides() const override;
        virtual std::vector<std::shared_ptr<Product>>& needs() override; // returns empty

        // TODO: undefine these empty ones, it's not compiling without them for some reason
        virtual validity is_valid_file(const std::string& file, FileReaderHints& hints){(void)file; (void)hints; return invalid;}
        validity is_valid_file(const std::string& file){
            FileReaderHints hints{file};
            return is_valid_file(file, hints);
        }

        virtual std::unique_ptr<FileReader> initialize_reader(DataProviderConfiguration& configuration, const boost::filesystem::path& path) = 0;

        virtual TileParameters read_next_tile(DataProviderConfiguration& configuration, DataRecord& record){(void)configuration; (void)record; return {0,0}; } // true on finished? For now.
};


} // namespace focs

#endif // FOCS_FILEREADER

