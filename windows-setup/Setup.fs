(*
 * Permission is granted to copy, distribute and/or modify the documents
 * in this directory and its subdirectories unless otherwise stated under
 * the terms of the GNU Free Documentation License, Version 1.1 or any later version 
 * published by the Free Software Foundation; with no Invariant Sections, 
 * no Front-Cover Texts and no Back-Cover Texts. A copy of the license 
 * is available at the website of the GNU Project.
 * The programs and code snippets in this directory and its subdirectories
 * are free software; you can redistribute them and/or modify it under the 
 * terms of the GNU General Public License as published by the Free Software 
 * Foundation; either version 2 of the License, or (at your option) any later
 * version. This code is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * Author Marco M. Mosca, email: marcomichele.mosca@gmail.com
*)
open System.IO
open System.Net.Http
open System.Diagnostics
open System

let download (output_dir:string) (http_client: HttpClient) (uri: string) = 
    let path = Path.Combine(output_dir, Path.GetFileName(uri))
    async {
        let fileBytes = http_client.GetByteArrayAsync(uri).Result
        File.WriteAllBytes(path, fileBytes);
    } |> Async.RunSynchronously
    path

let start_process (proc: string) (args: string array) =
    let startInfo = ProcessStartInfo()
    for a in args do
        startInfo.ArgumentList.Add(a)
    startInfo.FileName <- proc
    startInfo.RedirectStandardOutput <- true
    startInfo.RedirectStandardError <- true
    startInfo.UseShellExecute <- false
    startInfo.CreateNoWindow <- true
    let proc = new Process(StartInfo = startInfo)
    proc.Start() |> ignore
    (proc.StandardOutput.ReadToEnd()) |> ignore
    printfn "%s" (proc.StandardError.ReadToEnd())
    proc.WaitForExit()

let configure (src_dir:string) (variables: string array) = 
    let build_dir = Path.Combine(src_dir, "build")
    [| 
        "-B"; build_dir; 
        "-S"; src_dir;
        "-L";
        "-G"; "Visual Studio 17 2022";
        sprintf @"-DCMAKE_INSTALL_PREFIX=%s" (Path.Combine(src_dir, "install"));
        @"-DCMAKE_C_COMPILER=cl";
        @"-DBUILD_EXAMPLES=OFF";
        @"-DBUILD_TESTING=OFF";
    |]
    |> Array.append variables
    |> start_process "cmake"
     
    build_dir

let install (build_dir: string) =
    [|
        "--build"; build_dir;
        "--target"; "install";
        "--config"; "Release"
    |]
    |> start_process "cmake"

let unzip (output_folder: string) (file: string) =
    [| "x"; "-tzip"; file; sprintf "-o%s" output_folder |]
    |> start_process "7z"


[<EntryPoint>]
let main args =
    let ext_dir = Path.GetFullPath(@"..\External")
    assert Directory.Exists(ext_dir)
    let boost_dir = Path.Combine(ext_dir, @"boost_1_69_0")

    assert Directory.Exists(ext_dir)
    let src_nobuild_deps = [|
        (@"https://gitlab.com/libeigen/eigen/-/archive/3.3.7/eigen-3.3.7.zip", @"eigen-3.3.7");
        (@"https://github.com/project-gemmi/gemmi/archive/refs/tags/v0.3.3.zip", @"gemmi-0.3.3");
        (@"https://github.com/jlblancoc/nanoflann/archive/refs/tags/v1.3.1.zip", @"nanoflann-1.3.1");
        (@"https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/npackd/org.xmlsoft.LibXML-2.7.8.zip", @"libxml2-2.7.8.win32");
        
    |]

    let src_deps = [|
        (@"https://www.zlib.net/zlib1213.zip", @"zlib-1.2.13");
        (@"https://github.com/openbabel/openbabel/archive/refs/tags/openbabel-3-1-1.zip", @"openbabel-openbabel-3-1-1");
    |]
    let http = new HttpClient()
    for dep, folder_name in src_nobuild_deps do
        printfn "%s - %s" folder_name "Downloading..."
        let path_to_file = download ext_dir http dep
        printfn "%s - %s" folder_name "Unzipping..."
        unzip ext_dir path_to_file
        File.Delete(path_to_file)
    for dep, folder_name in src_deps do
        printfn "%s - %s" folder_name "Downloading..."
        let path_to_file = download ext_dir http dep
        printfn "%s - %s" folder_name "Unzipping..."
        unzip ext_dir path_to_file
        File.Delete(path_to_file)
        printfn "%s - %s" folder_name "Configuring..."
        let build_dir = configure (Path.Combine(ext_dir, folder_name)) [||]
        printfn "%s - %s" folder_name "Installing..."
        install build_dir


    let small_hpp = Path.Combine(ext_dir, "gemmi-0.3.3/include/gemmi/small.hpp")
    assert File.Exists(small_hpp)
    let new_lines = [|
        @"  if (len > 1) {";
        @"    if (label[len + 1] == '-')";
        @"      dest->charge = -dest->charge;";
        @"  }"
    |]
    let small_lines = 
        System.IO.File.ReadAllLines(small_hpp)
        |> Array.indexed 
        |> Array.filter (fun (i,x) -> i <> 69 && i <> 70)
        |> Array.map (fun (i,x) -> x) 
        |> Array.insertManyAt 69 new_lines
    File.Delete(small_hpp)
    File.WriteAllLines(small_hpp, small_lines)

    let babel_src = Path.Combine(ext_dir, @"openbabel-openbabel-3-1-1")
    let babel_dir = Path.Combine(ext_dir, @"openbabel-3-1-1\install")
    let babel_dir_include = Path.Combine(babel_dir, @"include")
    let babel_dir_bin = Path.Combine(babel_dir, @"bin")
    let babel_dir_data = Path.Combine(babel_dir, @"data")
    Directory.CreateDirectory(babel_dir) |> ignore
    Directory.Move(Path.Combine(babel_src, @"include"), babel_dir_include)
    File.Move(Path.Combine(babel_src, @"build\include\openbabel\babelconfig.h"), Path.Combine(babel_dir_include, @"openbabel\babelconfig.h"))
    Directory.Move(Path.Combine(babel_src, @"build\bin\Release"), babel_dir_bin)
    Directory.Move(Path.Combine(babel_src, @"build\src\Release\openbabel-3.lib"), Path.Combine(babel_dir_bin, @"openbabel-3.lib"))
    Directory.Move(Path.Combine(babel_src, @"build\src\Release\openbabel-3.exp"), Path.Combine(babel_dir_bin, @"openbabel-3.exp"))
    Directory.Move(Path.Combine(babel_src, @"data"), babel_dir_data)

    let zlib_dir = Path.Combine(ext_dir, @"zlib-1.2.13\install")
    let libxml2_dir = Path.Combine(ext_dir, @"libxml2-2.7.8.win32")
    for path_value in [|
        Path.Combine(boost_dir, @"lib64-msvc-14.1");
        Path.Combine(libxml2_dir, @"bin");
        Path.Combine(zlib_dir, @"bin"); 
        babel_dir_bin |] do
        assert Directory.Exists(path_value)
        let mutable path = Environment.GetEnvironmentVariable("Path", EnvironmentVariableTarget.User);
        if not (path.Contains(path_value)) then
            path <- path + sprintf ";%s" path_value;
            Environment.SetEnvironmentVariable("Path", path, EnvironmentVariableTarget.User);
    Environment.SetEnvironmentVariable("BABEL_DATADIR", babel_dir_data, EnvironmentVariableTarget.User);
    Environment.SetEnvironmentVariable("BABEL_LIBDIR", babel_dir_bin, EnvironmentVariableTarget.User);


    printfn "%s - %s" "Project AMD" "Configuring..."
    let amd_build_dir = configure (Path.GetFullPath(@"..\")) [||]
    printfn "%s - %s" "Project AMD" "Installing..."
    install amd_build_dir
    let mutable path = Environment.GetEnvironmentVariable("Path", EnvironmentVariableTarget.User);
    let amd_bin = sprintf @";%s\bin" (Path.GetFullPath(@"..\install"))
    if not (path.Contains(amd_bin)) then
        path <- path + amd_bin
        Environment.SetEnvironmentVariable("Path", path, EnvironmentVariableTarget.User);
    printfn "%s - %s" "Project AMD" "Completed!"
    printfn "%s - %s" "Project AMD" "Restart your command prompt to use the executables!"
    0